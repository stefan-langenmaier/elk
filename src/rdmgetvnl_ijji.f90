
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmgetvnl_ijji
! !INTERFACE:
subroutine rdmgetvnl_ijji(ikp,vnlijji)
! !USES:
use modrdm
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced k-point set (in,integer)
!   vnlijji : non-local Coulomb matrix elements (out,complex(nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Gets non-local Coulomb matrix elements of the type $(i-jj-i)$ from the file
!   {\tt RDMVNL\_IJJI.OUT}.
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(out) :: vnlijji(nstsv,nstsv,nkpt)
integer recl,iostat
! allocate arrays
inquire(iolength=recl) vnlijji
!$OMP CRITICAL
open(97,file='RDMVNL_IJJI.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl,iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(rdmgetvnl_ijji): error opening RDMVNL_IJJI.OUT")')
  write(*,*)
  stop
end if
read(97,rec=ikp) vnlijji
close(97)
!$OMP END CRITICAL
return
end subroutine
!EOC
