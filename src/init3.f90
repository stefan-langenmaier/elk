
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init3
use modmain
implicit none
! local variables
integer ig

!-----------------------!
!     RPA variables     !
!-----------------------!
! frequencies
nwrpa=1
if (allocated(wrpa)) deallocate(wrpa)
allocate(wrpa(nwrpa))
wrpa(1)=cmplx(0.d0,swidth,8)
! G-vectors
ngrpa=1
do ig=ngvec,1,-1
  if (gc(ig).lt.gmaxrpa) then
    ngrpa=ig
    goto 10
  end if
end do
10 continue

return
end subroutine

