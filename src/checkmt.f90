
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: checkmt
! !INTERFACE:
subroutine checkmt
! !USES:
use modmain
! !DESCRIPTION:
!   Checks for muffin-tins which are too close together. If any such muffin-tins
!   are found then the radii of their associated atomic species are adjusted so
!   that the minimum distance between their surfaces is {\tt rmtdelta}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!   Modified, October 2011 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,js
real(8) dmin,t1,t2
! automatic arrays
real(8) rmt0(nspecies)
rmt0(1:nspecies)=rmt(1:nspecies)
10 continue
! find the minimum distance between muffin-tin surfaces
call mtdmin(is,js,dmin)
! adjust muffin-tin radii if required
if (dmin.lt.rmtdelta-epslat) then
  t1=rmt(is)+rmt(js)
  t2=(t1+dmin-rmtdelta)/t1
  rmt(is)=rmt(is)*t2
  if (is.ne.js) rmt(js)=rmt(js)*t2
  goto 10
end if
do is=1,nspecies
  if (rmt(is).lt.0.25d0) then
    write(*,*)
    write(*,'("Error(checkmt): muffin-tin radius too small for species ",I4,&
     &" (",A,")")') is,trim(spsymb(is))
    write(*,'(" Radius : ",G18.10)') rmt(is)
    write(*,*)
    stop
  end if
  if (rmt(is).lt.rmt0(is)) then
    write(*,*)
    write(*,'("Info(checkmt): muffin-tin radius reduced from ",G18.10," to ",&
     &G18.10)') rmt0(is),rmt(is)
    write(*,'(" for species ",I4," (",A,")")') is,trim(spsymb(is))
  end if
end do
return
end subroutine
!EOC
