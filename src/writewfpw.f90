
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writewfpw
use modmain
use modmpi
implicit none
! local variables
logical twfpwh
integer ik,recl
! allocatable arrays
complex(8), allocatable :: wfpw(:,:,:)
complex(8), allocatable :: wfpwh(:,:,:,:,:)
! initialise global variables
call init0
call init1
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! require high planewave wavefunctions or not
if (task.eq.135) then
  twfpwh=.false.
else
  twfpwh=.true.
end if
! delete existing WFPW.OUT and WFPWH.OUT
if (mp_mpi) then
  open(50,file='WFPW.OUT')
  close(50,status='DELETE')
  if (twfpwh) then
    open(51,file='WFPWH.OUT')
    close(51,status='DELETE')
  end if
end if
! synchronise MPI processes
call mpi_barrier(mpi_comm_world,ierror)
! determine the record length and open WFPW.OUT
allocate(wfpw(ngkmax,nspinor,nstsv))
inquire(iolength=recl) vkl(:,1),ngkmax,nspinor,nstsv,wfpw
deallocate(wfpw)
open(50,file='WFPW.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
! determine the record length and open WFPWH.OUT
if (twfpwh) then
  allocate(wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  inquire(iolength=recl) vkl(:,1),lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv,wfpwh
  deallocate(wfpwh)
  open(51,file='WFPWH.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
   recl=recl)
end if
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfpw,wfpwh)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(wfpw(ngkmax,nspinor,nstsv))
  if (twfpwh) allocate(wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
!$OMP CRITICAL
  write(*,'("Info(writewfpw): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! generate the plane wave wavefunctions
  call genwfpw(twfpwh,vkl(:,ik),ngk(1,ik),igkig(:,1,ik),vgkl(:,:,1,ik), &
   gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),wfpw,wfpwh)
!$OMP CRITICAL
  write(50,rec=ik) vkl(:,ik),ngkmax,nspinor,nstsv,wfpw
  if (twfpwh) write(51,rec=ik) vkl(:,ik),lmmaxvr,nrcmtmax,natmtot,nspinor, &
   nstsv,wfpwh
!$OMP END CRITICAL
  deallocate(wfpw)
  if (twfpwh) deallocate(wfpwh)
end do
!$OMP END DO
!$OMP END PARALLEL
close(50)
if (twfpwh) close(51)
! synchronise MPI processes
call mpi_barrier(mpi_comm_world,ierror)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writewfpw):")')
  write(*,'(" low plane wave wavefunctions written to WFPW.OUT")')
  if (twfpwh) write(*,'(" high plane wave wavefunctions written to WFPWH.OUT")')
end if
return
end subroutine

