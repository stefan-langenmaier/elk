
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getvnlijji
! !INTERFACE:
subroutine getvnlijji(ikp,vnlijji)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced k-point set (in,integer)
!   vnlijji : non-local Coulomb matrix elements (out,real(nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Gets non-local Coulomb matrix elements of the type $(i-jj-i)$ from the file
!   {\tt VNLIJJI.OUT}.
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(out) :: vnlijji(nstsv,nstsv,nkpt)
! local variables
integer recl,iostat
! determine record length
inquire(iolength=recl) vnlijji
!$OMP CRITICAL
open(95,file='VNLIJJI.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl,iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(getvnlijji): error opening file VNLIJJI.OUT")')
  write(*,*)
  stop
end if
read(95,rec=ikp) vnlijji
close(95)
!$OMP END CRITICAL
return
end subroutine
!EOC

