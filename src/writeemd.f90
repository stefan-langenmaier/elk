
! Copyright (C) 2012 S. Dugdale, D. Ernsting and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeemd
use modmain
use modmpi
implicit none
! local variables
logical reduceh0
integer ik,ih,recl
integer ist,ispn
real(8) sum,t1
complex(8) zt1
! allocatable arrays
real(8), allocatable :: vhkc(:,:)
real(8), allocatable :: emd(:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfhk(:,:,:)
if (spinsprl) then
  write(*,*)
  write(*,'("Error(writeemd): electron momentum density not available for &
   &spin-spirals")')
  write(*,*)
  stop
end if
! initialise universal variables
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
! get the occupancies from file
do ik=1,nkpt
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! generate the H-vectors without reduction
reduceh0=reduceh
reduceh=.false.
call genhvec
reduceh=reduceh0
! delete existing EMD.OUT
if (mp_mpi) then
  open(85,file='EMD.OUT')
  close(85,status='DELETE')
end if
! synchronise MPI processes
call mpi_barrier(mpi_comm_world,ierror)
allocate(emd(nhvec))
inquire(iolength=recl) vkl(:,1),nhvec,emd
deallocate(emd)
open(85,file='EMD.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
! loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vhkc,emd,wfmt,wfir,wfhk) &
!$OMP PRIVATE(ih,sum,ist,ispn,zt1,t1)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(vhkc(3,nhvec))
  allocate(emd(nhvec))
  allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir(ngkmax,nspinor,nstsv))
  allocate(wfhk(nhvec,nspinor,nstsv))
!$OMP CRITICAL
  write(*,'("Info(writeemd): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! generate the second-variational wavefunctions
  call genwfsvp(.true.,.true.,vkl(:,ik),wfmt,ngkmax,wfir)
! generate the H+k-vectors
  do ih=1,nhvec
    vhkc(:,ih)=vhc(:,ih)+vkc(:,ik)
  end do
! Fourier transform the wavefunctions
  call genwfhk(ik,vhkc,wfmt,wfir,wfhk)
! loop over all H+k-vectors
  do ih=1,nhvec
! sum over occupied states and spins
    sum=0.d0
    do ist=1,nstsv
      do ispn=1,nspinor
        zt1=wfhk(ih,ispn,ist)
        t1=dble(zt1)**2+aimag(zt1)**2
        sum=sum+occsv(ist,ik)*t1
      end do
    end do
    emd(ih)=sum
  end do
!$OMP CRITICAL
  write(85,rec=ik) vkl(:,ik),nhvec,emd
!$OMP END CRITICAL
  deallocate(vhkc,emd,wfmt,wfir,wfhk)
end do
!$OMP END DO
!$OMP END PARALLEL
close(85)
! synchronise MPI processes
call mpi_barrier(mpi_comm_world,ierror)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writeemd): electron momentum density written to EMD.OUT")')
  write(*,'(" for all H+k-vectors up to |H| < hmax")')
end if
return
end subroutine

