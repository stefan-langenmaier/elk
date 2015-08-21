
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potnucl(ptnucl,nr,r,zn,vn)
implicit none
! arguments
logical, intent(in) :: ptnucl
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: zn
real(8), intent(out) :: vn(nr)
! local variables
integer ir
real(8) rn,t1,t2
! external functions
real(8) radnucl
external radnucl
if (zn.eq.0.d0) then
  vn(:)=0.d0
  return
end if
if (ptnucl) then
! nucleus is taken to be a point particle
  vn(:)=zn/r(:)
else
! approximate nuclear radius
  rn=radnucl(zn)
  t1=zn/(2.d0*rn**3)
  t2=3.d0*rn**2
  do ir=1,nr
    if (r(ir).lt.rn) then
      vn(ir)=t1*(t2-r(ir)**2)
    else
      vn(ir)=zn/r(ir)
    end if
  end do
end if
return
end subroutine

