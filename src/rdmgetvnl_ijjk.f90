
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmgetvnl_ijjk
! !INTERFACE:
subroutine rdmgetvnl_ijjk(ikp,vnlijjk)
! !USES:
use modrdm
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced k-point set (in,integer)
!   vnlijjk : non-local Coulomb matrix elements
!             (out,complex(nstsv,nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Gets Coulomb matrix elements of the type $(i-jj-k)$ from the file
!   {\tt RDMVNL\_IJJK.OUT}.
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vnlijjk(nstsv,nstsv,nstsv,nkpt)
integer recl,iostat
inquire(iolength=recl) vnlijjk
!$OMP CRITICAL
open(98,file='RDMVNL_IJJK.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl,iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(rdmgetvnl_ijjk): error opening RDMVNL_IJJK.OUT")')
  write(*,*)
  stop
end if
read(98,rec=ikp) vnlijjk
close(98)
!$OMP END CRITICAL
return
end subroutine
!EOC
