
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmwputvnl_ijji
! !INTERFACE:
subroutine rdmputvnl_ijji
! !USES:
use modrdm
use modmain
! !DESCRIPTION:
!   Generates non-local Coulomb matrix elements of the type $(i-jj-i)$ and
!   outputs them to the file {\tt RDMVNL\_IJJI.OUT}.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! allocatable arrays
real(8), allocatable :: vnlijji(:,:,:)
integer recl,ik
! allocate arrays
allocate(vnlijji(nstsv,nstsv,nkpt))
inquire(iolength=recl) vnlijji
deallocate(vnlijji)
open(94,file='RDMVNL_IJJI.OUT',action='WRITE',form='UNFORMATTED', &
 access='DIRECT',status='REPLACE',recl=recl)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE (vnlijji)
!$OMP DO
do ik=1,nkptnr
  allocate(vnlijji(nstsv,nstsv,nkpt))
!$OMP CRITICAL
  write(*,'("Info(rdmputvnl_IJJI): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
  call rdmgenvnl_ijji(ik,vnlijji)
!$OMP CRITICAL
  write(94,rec=ik) vnlijji
!$OMP END CRITICAL
  deallocate(vnlijji)
end do
!$OMP END DO
!$OMP END PARALLEL
close(94)
return
end subroutine
!EOC
