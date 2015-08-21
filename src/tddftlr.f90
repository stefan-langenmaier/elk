
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddftlr
use modmain
use modmpi
implicit none
! local variables
integer ik,ig,jg,iw
integer iq,isym,it,n
integer info1,info2
real(8) vecqc(3),vgqc(3)
real(8) v(3),t1
complex(8) vfxcp,zt1
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: gqc(:)
complex(8), allocatable :: expqmt(:,:,:)
complex(8), allocatable :: vchi0(:,:,:),vfxc(:,:,:)
complex(8), allocatable :: eps0(:,:,:),eps(:,:,:)
complex(8), allocatable :: a(:,:),work(:)
! initialise global variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
! read Fermi energy from a file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! allocate local arrays
allocate(ipiv(ngrpa))
allocate(gqc(ngrpa))
allocate(expqmt(lmmaxvr,nrcmtmax,natmtot))
allocate(vchi0(ngrpa,ngrpa,nwrpa),vfxc(ngrpa,ngrpa,nwrpa))
allocate(eps0(ngrpa,ngrpa,nwrpa),eps(ngrpa,ngrpa,nwrpa))
allocate(a(ngrpa,ngrpa),work(ngrpa))
! generate the exp(iG.r) functions for all the MBPT G-vectors
call genexpigr
! check q-vector is commensurate with k-point grid
v(:)=dble(ngridk(:))*vecql(:)
v(:)=abs(v(:)-nint(v(:)))
if ((v(1).gt.epslat).or.(v(2).gt.epslat).or.(v(3).gt.epslat)) then
  write(*,*)
  write(*,'("Error(tddftlr): q-vector incommensurate with k-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" vecql : ",3G18.10)') vecql
  write(*,*)
  stop
end if
call r3mv(bvec,vecql,vecqc)
call findqpt(vecql,isym,iq)
! generate the G+q vector lengths
do ig=1,ngrpa
  vgqc(:)=vgc(:,ig)+vecqc(:)
  gqc(ig)=sqrt(vgqc(1)**2+vgqc(2)**2+vgqc(3)**2)
end do
! generate the exp(iq.r) function
call genexpmt(vecqc(:),expqmt)
! compute v^1/2 chi0 v^1/2 (the symmetric version of v chi0)
vchi0(:,:,:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL
  write(*,'("Info(tddftlr): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
! compute v^1/2 chi0 v^1/2
  call genvchi0(ik,optcomp(1,1),scissor,vecql,gqc,expqmt,vchi0)
end do
!$OMP END DO
!$OMP END PARALLEL
! add epsinv from each process and redistribute
if (np_mpi.gt.1) then
  n=ngrpa*ngrpa*nwrpa
  call mpi_allreduce(mpi_in_place,vchi0,n,mpi_double_complex,mpi_sum, &
   mpi_comm_world,ierror)
end if
! calculate symmetric epsilon = 1 - v^1/2 chi0 v^1/2
do ig=1,ngrpa
  do jg=1,ngrpa
    eps0(ig,jg,:)=-vchi0(ig,jg,:)
    eps(ig,jg,:)=vchi0(ig,jg,:)
  end do
  eps0(ig,ig,:)=eps0(ig,ig,:)+1.d0
  eps(ig,ig,:)=eps(ig,ig,:)+1.d0
end do
vfxcp=0.d0
it=0
10 continue
! compute vchi0 v^(-1/2) f_xc v^(-1/2) vchi0
call genvfxc(gqc,vchi0,eps0,eps,vfxc)
! begin loop over frequencies
do iw=1,nwrpa
! compute 1 - v^1/2 chi0 v^1/2 - v^(-1/2) f_xc v^(-1/2) vchi0
  a(:,:)=eps0(:,:,iw)-vfxc(:,:,iw)
! invert this matrix
  call zgetrf(ngrpa,ngrpa,a,ngrpa,ipiv,info1)
  call zgetri(ngrpa,a,ngrpa,ipiv,work,ngrpa,info2)
  if ((info1.ne.0).or.(info2.ne.0)) then
    write(*,*)
    write(*,'("Error(tddftlr): unable to invert epsilon")')
    write(*,'(" for MBPT frequency ",I6)') iw
    write(*,*)
    stop
  end if
! left multiply by v^1/2 chi0 v^1/2
  call zgemm('N','N',ngrpa,ngrpa,ngrpa,zone,vchi0(:,:,iw),ngrpa,a,ngrpa,zzero, &
   eps(:,:,iw),ngrpa)
! compute epsilon = 1 + v^1/2 chi v^1/2
  do ig=1,ngrpa
    eps(ig,ig,iw)=1.d0+eps(ig,ig,iw)
  end do
end do
! bootstrap f_xc
if (fxctype(1).eq.210) then
  it=it+1
  if (it.gt.500) then
    write(*,*)
    write(*,'("Error(tddftlr): bootstrap kernel failed to converge")')
    write(*,*)
    stop
  end if
! check for convergence
  t1=abs(vfxcp)-abs(vfxc(1,1,1))
  vfxcp=vfxc(1,1,1)
  if (abs(t1).gt.1.d-8) goto 10
end if
! write G = G' = 0 components to file
open(50,file="EPSILON_TDDFT.OUT",action='WRITE',form='FORMATTED')
open(51,file="EELS_TDDFT.OUT",action='WRITE',form='FORMATTED')
do iw=1,nwrpa
  zt1=1.d0/eps(1,1,iw)
  write(50,'(3G18.10)') dble(wrpa(iw)),dble(zt1)
  write(51,'(3G18.10)') dble(wrpa(iw)),-dble(eps(1,1,iw))
end do
write(50,*)
write(51,*)
do iw=1,nwrpa
  zt1=1.d0/eps(1,1,iw)
  write(50,'(3G18.10)') dble(wrpa(iw)),aimag(zt1)
  write(51,'(3G18.10)') dble(wrpa(iw)),-aimag(eps(1,1,iw))
end do
close(50)
close(51)
write(*,*)
write(*,'("Info(tddftlr):")')
write(*,'(" dielectric tensor written to EPSILON_TDDFT.OUT")')
write(*,'(" electron loss function written to EELS_TDDFT.OUT")')
write(*,'(" for component i, j = ",I1)') optcomp(1,1)
write(*,'(" q-vector (lattice coordinates) : ")')
write(*,'(3G18.10)') vecql
write(*,'(" q-vector length : ",G18.10)') gqc(1)
deallocate(ipiv,gqc,expqmt,vchi0,vfxc)
deallocate(eps0,eps,a,work)
! deallocate global exp(iG.r) arrays
deallocate(expgmt,expgir)
return
end subroutine
