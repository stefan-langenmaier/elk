
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddftlr
use modmain
use modmpi
implicit none
! local variables
integer ik,ig,jg,iw,n
integer info1,info2
complex(8) zt1
! allocatable arrays
integer, allocatable :: ipiv(:)
complex(8), allocatable :: expqmt(:,:,:)
complex(8), allocatable :: vchi0(:,:,:),fxc(:,:,:)
complex(8), allocatable :: eps0(:,:,:),eps(:,:,:)
complex(8), allocatable :: vc(:,:),e(:,:),vce(:,:)
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
allocate(expqmt(lmmaxvr,nrcmtmax,natmtot))
allocate(vchi0(nwrpa,ngrpa,ngrpa),fxc(nwrpa,ngrpa,ngrpa))
allocate(eps0(nwrpa,ngrpa,ngrpa),eps(nwrpa,ngrpa,ngrpa))
allocate(vc(ngrpa,ngrpa),e(ngrpa,ngrpa),vce(ngrpa,ngrpa))
allocate(a(ngrpa,ngrpa),work(ngrpa))
! generate the exp(iG.r) functions for all the RPA G-vectors
call genexpigr
! generate the exp(iq.r) function for q = 0
expqmt(:,:,:)=1.d0
! compute v^1/2 chi0 v^1/2 for q = 0 (the symmetric version of v chi0)
vchi0(:,:,:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call genvchi0(iq0,ik,gc,expqmt,vchi0)
end do
!$OMP END DO
!$OMP END PARALLEL
! add epsinv from each process and redistribute
if (np_mpi.gt.1) then
  n=nwrpa*ngrpa*ngrpa
  call mpi_allreduce(mpi_in_place,vchi0,n,mpi_double_complex,mpi_sum, &
   mpi_comm_world,ierror)
end if
! calculate symmetric epsilon = 1 - v^1/2 chi0 v^1/2
do ig=1,ngrpa
  do jg=1,ngrpa
    eps0(:,ig,jg)=-vchi0(:,ig,jg)
  end do
  eps0(:,ig,ig)=eps0(:,ig,ig)+1.d0
end do
! compute vchi0 v^1/2 fxc v^1/2 vchi0
call genfxc(vchi0,fxc)
! begin loop over frequencies
do iw=1,nwrpa
  vc(:,:)=vchi0(iw,:,:)
  e(:,:)=eps0(iw,:,:)
! left multiply eps0 by vchi0 to get v^1/2 chi0 v^1/2 - v^1/2 chi0 v chi0 v^1/2
  call zgemm('N','N',ngrpa,ngrpa,ngrpa,zone,vc,ngrpa,e,ngrpa,zzero,vce,ngrpa)
! subtract v^1/2 chi0 fxc chi0 v^1/2
  vce(:,:)=vce(:,:)-fxc(iw,:,:)
! invert this matrix
  call zgetrf(ngrpa,ngrpa,vce,ngrpa,ipiv,info1)
  call zgetri(ngrpa,vce,ngrpa,ipiv,work,ngrpa,info2)
  if ((info1.ne.0).or.(info2.ne.0)) then
    write(*,*)
    write(*,'("Error(tddftlr): unable to invert epsilon")')
    write(*,'(" for RPA frequency ",I6)') iw
    write(*,*)
    stop
  end if
! compute v^1/2 chi v^1/2 = vchi0 vce vchi0
  call zgemm('N','N',ngrpa,ngrpa,ngrpa,zone,vc,ngrpa,vce,ngrpa,zzero,a,ngrpa)
  call zgemm('N','N',ngrpa,ngrpa,ngrpa,zone,a,ngrpa,vc,ngrpa,zzero,e,ngrpa)
! compute epsilon = 1 + v^1/2 chi v^1/2
  eps(iw,:,:)=e(:,:)
  do ig=1,ngrpa
    eps(iw,ig,ig)=1.d0+eps(iw,ig,ig)
  end do
end do
! write G = G' = 0 components to file
open(50,file="EPSILON_TDDFT.OUT",action='WRITE',form='FORMATTED')
do iw=1,nwrpa
  zt1=1.d0/eps(iw,1,1)
  write(50,'(3G18.10)') dble(wrpa(iw)),dble(zt1)
end do
write(50,*)
do iw=1,nwrpa
  zt1=1.d0/eps(iw,1,1)
  write(50,'(3G18.10)') dble(wrpa(iw)),aimag(zt1)
end do
close(50)
write(*,*)
write(*,'("Info(tddftlr):")')
write(*,'(" 1/3 trace of dielectric tensor written to EPSILON_TDDFT.OUT")')
write(*,*)
deallocate(ipiv,expqmt,vchi0,fxc)
deallocate(eps0,eps,vc,e,vce,a,work)
! deallocate global exp(iG.r) arrays
deallocate(expgmt,expgir)
return
end subroutine

