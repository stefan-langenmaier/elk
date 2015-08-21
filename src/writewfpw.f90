
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writewfpw
use modmain
use modmpi
implicit none
! local variables
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
! delete existing WFPW.OUT and WFPWH.OUT
if (mp_mpi) then
  open(50,file='WFPW.OUT')
  close(50,status='DELETE')
  open(51,file='WFPWH.OUT')
  close(51,status='DELETE')
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
allocate(wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
inquire(iolength=recl) vkl(:,1),lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv,wfpwh
deallocate(wfpwh)
open(51,file='WFPWH.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfpw,wfpwh)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(wfpw(ngkmax,nspinor,nstsv))
  allocate(wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
!$OMP CRITICAL
  write(*,'("Info(writewfpw): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! generate the plane wave wavefunctions
  call genwfpw(vkl(:,ik),ngk(1,ik),igkig(:,1,ik),vgkl(:,:,1,ik),gkc(:,1,ik), &
   tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),wfpw,wfpwh)
!$OMP CRITICAL
  write(50,rec=ik) vkl(:,ik),ngkmax,nspinor,nstsv,wfpw
  write(51,rec=ik) vkl(:,ik),lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv,wfpwh
!$OMP END CRITICAL
  deallocate(wfpw,wfpwh)
end do
!$OMP END DO
!$OMP END PARALLEL
close(50)
close(51)
! synchronise MPI processes
call mpi_barrier(mpi_comm_world,ierror)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writewfpw): low and high plane wave wavefunctions written &
   &to")')
  write(*,'(" WFPW.OUT and WFPWH.OUT, respectively")')
  write(*,*)
end if
return
end subroutine

