
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genffacg
use modmain
implicit none
! local variables
integer is,ig
real(8) t1,t2
if (allocated(ffacg)) deallocate(ffacg)
allocate(ffacg(ngrtot,nspecies))
t1=fourpi/omega
do is=1,nspecies
  ffacg(1,is)=(t1/3.d0)*rmt(is)**3
  do ig=2,ngrtot
    t2=gc(ig)*rmt(is)
    ffacg(ig,is)=t1*(sin(t2)-t2*cos(t2))/(gc(ig)**3)
  end do
end do
return
end subroutine

