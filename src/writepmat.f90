
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmat
! !INTERFACE:
subroutine writepmat
! !USES:
use modmain
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
integer ik,ispn,recl
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
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
! find the record length
allocate(pmat(3,nstsv,nstsv))
inquire(iolength=recl) pmat
deallocate(pmat)
open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 status='REPLACE',recl=recl)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv,pmat,ispn)
!$OMP DO
do ik=1,nkpt
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(evecfv(nmatmax,nstfv))
  allocate(evecsv(nstsv,nstsv))
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
! calculate the momentum matrix elements
  call genpmat(ngk(:,ik),igkig(:,:,ik),vgkc(:,:,:,ik),apwalm,evecfv,evecsv,pmat)
! write the matrix elements to direct-access file
!$OMP CRITICAL
  write(50,rec=ik) pmat
!$OMP END CRITICAL
  deallocate(apwalm,evecfv,evecsv,pmat)
end do
!$OMP END DO
!$OMP END PARALLEL
close(50)
write(*,*)
write(*,'("Info(writepmat):")')
write(*,'(" momentum matrix elements written to file PMAT.OUT")')
write(*,*)
end subroutine
!EOC

