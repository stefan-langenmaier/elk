
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine epsinv_rpa
use modmain
use modmpi
implicit none
! local variables
integer iq,ik,ig,jg,iw,n
integer info1,info2,recl
real(8) vgqc(3)
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: gqc(:)
complex(8), allocatable :: expqmt(:,:,:)
complex(8), allocatable :: epsinv(:,:,:)
complex(8), allocatable :: work(:)
! initialise global variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! generate the exp(iG.r) functions for all the MBPT G-vectors
call genexpigr
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
allocate(gqc(ngrpa))
allocate(expqmt(lmmaxvr,nrcmtmax,natmtot))
allocate(epsinv(ngrpa,ngrpa,nwrpa))
if (mp_mpi) then
! determine the record length for EPSINV_RPA.OUT
  inquire(iolength=recl) vql(:,1),ngrpa,nwrpa,epsinv
! open EPSINV_RPA.OUT
  open(50,file='EPSINV_RPA.OUT',action='WRITE',form='UNFORMATTED', &
   access='DIRECT',status='REPLACE',recl=recl)
end if
! loop over q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(epsinv_rpa): ",I6," of ",I6," q-points")') iq,nqpt
! generate the G+q-vector lengths
  do ig=1,ngrpa
    vgqc(:)=vgc(:,ig)+vqc(:,iq)
    gqc(ig)=sqrt(vgqc(1)**2+vgqc(2)**2+vgqc(3)**2)
  end do
! generate the function exp(iq.r) in the muffin-tins
  call genexpmt(vqc(:,iq),expqmt)
! zero the dielectric function array
  epsinv(:,:,:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! compute v^1/2 chi0 v^1/2 and store in array epsinv
    call genvchi0(iq,ik,0,gqc,expqmt,epsinv)
  end do
!$OMP END DO
!$OMP END PARALLEL
! add epsinv from each process and redistribute
  if (np_mpi.gt.1) then
    n=ngrpa*ngrpa*nwrpa
    call mpi_allreduce(mpi_in_place,epsinv,n,mpi_double_complex,mpi_sum, &
     mpi_comm_world,ierror)
  end if
! negate and add delta(G,G')
  do ig=1,ngrpa
    do jg=1,ngrpa
      epsinv(ig,jg,:)=-epsinv(ig,jg,:)
    end do
    epsinv(ig,ig,:)=epsinv(ig,ig,:)+1.d0
  end do
!-------------------------------------!
!     invert epsilon over G-space     !
!-------------------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ipiv,work,info1,info2)
!$OMP DO
  do iw=1,nwrpa
    allocate(ipiv(ngrpa),work(ngrpa))
    call zgetrf(ngrpa,ngrpa,epsinv(:,:,iw),ngrpa,ipiv,info1)
    call zgetri(ngrpa,epsinv(:,:,iw),ngrpa,ipiv,work,ngrpa,info2)
    if ((info1.ne.0).or.(info2.ne.0)) then
      write(*,*)
      write(*,'("Error(epsinv_rpa): unable to invert epsilon")')
      write(*,'(" for q-point ",I6)') iq
      write(*,'(" and MBPT frequency ",I6)') iw
      write(*,*)
      stop
    end if
    deallocate(ipiv,work)
  end do
!$OMP END DO
!$OMP END PARALLEL
! write inverse RPA epsilon to EPSINV_RPA.OUT
  if (mp_mpi) write(50,rec=iq) vql(:,iq),ngrpa,nwrpa,epsinv
! end loop over q-points
end do
if (mp_mpi) close(50)
deallocate(gqc,expqmt,epsinv)
! deallocate global exp(iG.r) arrays
deallocate(expgmt,expgir)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(epsinv_rpa):")')
  write(*,'(" inverse RPA dielectric function, eps^(-1)(G,G'',q,w), written to &
   &EPSINV_RPA.OUT")')
end if
return
end subroutine

