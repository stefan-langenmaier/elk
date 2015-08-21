
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init3
use modmain
implicit none
! local variables
integer ig,iw
real(8) w1,w2,t1,t2

!--------------------------------------------------------!
!     many-body perturbation theory (MBPT) variables     !
!--------------------------------------------------------!
! G-vectors for 2-point correlators
ngrpa=1
do ig=ngvec,1,-1
  if (gc(ig).lt.gmaxrpa) then
    ngrpa=ig
    exit
  end if
end do
! frequencies for bosonic 2-point correlators
nwrpa=1
if (allocated(wrpa)) deallocate(wrpa)
if (task.eq.188) then
  nwrpa=nwdos
  allocate(wrpa(nwrpa))
  w1=max(wdos(1),0.d0)
  w2=max(wdos(2),w1)
  t1=(w2-w1)/dble(nwdos)
  do iw=1,nwdos
    t2=w1+t1*dble(iw-1)
    wrpa(iw)=cmplx(t2,swidth,8)
  end do
else
  nwrpa=1
  allocate(wrpa(nwrpa))
  wrpa(1)=cmplx(0.d0,swidth,8)
end if

return
end subroutine

