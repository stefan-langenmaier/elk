
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genffacg(is,ngv,ffacg)
use modmain
implicit none
! arguments
integer, intent(in) :: is
integer, intent(in) :: ngv
real(8), intent(out) :: ffacg(ngv)
! local variables
integer ig
real(8) t1,t2
t1=fourpi/omega
ffacg(1)=(t1/3.d0)*rmt(is)**3
do ig=2,ngv
  t2=gc(ig)*rmt(is)
  ffacg(ig)=t1*(sin(t2)-t2*cos(t2))/(gc(ig)**3)
end do
return
end subroutine

