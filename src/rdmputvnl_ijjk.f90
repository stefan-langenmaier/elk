
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmputvnl_ijjk
! !INTERFACE:
subroutine rdmputvnl_ijjk
! !USES:
use modrdm
use modmain
! !DESCRIPTION:
!   Generates non-local Coulomb matrix elements of the type $(i-jj-k)$ and
!   $(i-jj-i)$, and outputs them to the files {\tt RDMVNL\_IJJK.OUT} and
!   {\tt RDMVNL\_IJJI.OUT}, respectively.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! allocatable arrays
complex(8), allocatable :: vnlijjk(:,:,:,:)
real(8), allocatable :: vnlijji(:,:,:)
integer recl1,recl2,ik
allocate(vnlijjk(nstsv,nstsv,nstsv,nkpt))
allocate(vnlijji(nstsv,nstsv,nkpt))
inquire(iolength=recl1) vnlijjk
inquire(iolength=recl2) vnlijji
deallocate(vnlijjk,vnlijji)
open(95,file='RDMVNL_IJJK.OUT',action='WRITE',form='UNFORMATTED', &
 access='DIRECT',status='REPLACE',recl=recl1)
open(96,file='RDMVNL_IJJI.OUT',action='WRITE',form='UNFORMATTED', &
 access='DIRECT',status='REPLACE',recl=recl2)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vnlijjk,vnlijji)
!$OMP DO
do ik=1,nkptnr
  allocate(vnlijjk(nstsv,nstsv,nstsv,nkpt))
  allocate(vnlijji(nstsv,nstsv,nkpt))
!$OMP CRITICAL
  write(*,'("Info(rdmputvnl_ijjk): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
! calculate non-local matrix elements of the type (l-jj-k)
  call rdmgenvnl_ijjk(ik,vnlijji,vnlijjk)
!$OMP CRITICAL
  write(95,rec=ik) vnlijjk
  write(96,rec=ik) vnlijji
!$OMP END CRITICAL
  deallocate(vnlijjk,vnlijji)
end do
!$OMP END DO
!$OMP END PARALLEL
close(95)
close(96)
return
end subroutine
!EOC
