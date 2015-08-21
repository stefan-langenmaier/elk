
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmat
! !INTERFACE:
subroutine writepmat
! !USES:
use modmain
use modmpi
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT.OUT}.
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,ispn,recl
! allocatable arrays
integer, allocatable :: idx(:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: pmat(:,:,:)
! initialise universal variables
call init0
call init1
! read in the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! index to all states
allocate(idx(nstsv))
do ist=1,nstsv
  idx(ist)=ist
end do
! delete existing PMAT.OUT
if (mp_mpi) then
  open(85,file='PMAT.OUT')
  close(85,status='DELETE')
end if
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
! determine the record length for PMAT.OUT
allocate(pmat(3,nstsv,nstsv))
inquire(iolength=recl) vkl(:,1),nstsv,pmat
deallocate(pmat)
open(85,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecfv,evecsv,apwalm,wfmt,wfir,pmat,ispn)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir(ngkmax,nspinor,nstsv))
  allocate(pmat(3,nstsv,nstsv))
!$OMP CRITICAL
  write(*,'("Info(writepmat): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! get the eigenvectors from file
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! calculate the wavefunctions for all states
  call genwfsv(.true.,.true.,nstsv,idx,ngk(:,ik),igkig(:,:,ik),apwalm,evecfv, &
   evecsv,wfmt,ngkmax,wfir)
! calculate the momentum matrix elements
  call genpmat(ngk(:,ik),igkig(:,:,ik),vgkc(:,:,:,ik),wfmt,wfir,pmat)
! write the matrix elements to direct-access file
!$OMP CRITICAL
  write(85,rec=ik) vkl(:,ik),nstsv,pmat
!$OMP END CRITICAL
  deallocate(evecfv,evecsv,apwalm,wfmt,wfir,pmat)
end do
!$OMP END DO
!$OMP END PARALLEL
close(85)
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writepmat):")')
  write(*,'(" momentum matrix elements written to file PMAT.OUT")')
end if
deallocate(idx)
return
end subroutine
!EOC

