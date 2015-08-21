
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfmtinp
! !INTERFACE:
complex(8) function zfmtinp(tsh,nr,nri,r,zfmt1,zfmt2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh   : .true. if the functions are in spherical harmonics (in,logical)
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of points on the inner part of the muffin-tin (in,integer)
!   r     : radial mesh (in,real(nr))
!   zfmt1 : first complex muffin-tin function in spherical harmonics or
!           coordinates (in,complex(lmmaxvr,nr))
!   zfmt2 : second complex muffin-tin function in spherical harmonics or
!           coordinates (in,complex(lmmaxvr,nr))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions in the muffin-tin. In
!   other words, given two complex functions of the form
!   $$ f({\bf r})=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}(r)Y_{lm}
!    (\hat{\bf r}), $$
!   the function returns
!   $$ I=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}\int f_{lm}^{1*}(r)
!    f_{lm}^2(r)r^2\,dr\;. $$
!   Note that if {\tt tsh} is {\tt .false.} the functions are in spherical
!   coordinates rather than spherical harmonics. In this case $I$ is multiplied
!   by $4\pi/(l_{\rm max}+1)^2$.
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!   Modified, September 2013 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr)
complex(8), intent(in) :: zfmt1(lmmaxvr,nr),zfmt2(lmmaxvr,nr)
! local variables
integer nr0,ir0,ir
real(8) t1
complex(8) z1
! automatic arrays
real(8) fr1(nr),fr2(nr),gr(nr)
! external functions
complex(8) zdotc
external zdotc
if (tsh) then
! functions are in spherical harmonics
  do ir=1,nri
    z1=zdotc(lmmaxinr,zfmt1(:,ir),1,zfmt2(:,ir),1)*(r(ir)**2)
    fr1(ir)=dble(z1)
    fr2(ir)=aimag(z1)
  end do
  do ir=nri+1,nr
    z1=zdotc(lmmaxvr,zfmt1(:,ir),1,zfmt2(:,ir),1)*(r(ir)**2)
    fr1(ir)=dble(z1)
    fr2(ir)=aimag(z1)
  end do
  call fderiv(-1,nr,r,fr1,gr)
  t1=gr(nr)
  call fderiv(-1,nr,r,fr2,gr)
  zfmtinp=cmplx(t1,gr(nr),8)
else
! functions are in spherical coordinates
  if (nri.gt.0) then
    do ir=1,nri
      z1=zdotc(lmmaxinr,zfmt1(:,ir),1,zfmt2(:,ir),1)*(r(ir)**2)
      fr1(ir)=dble(z1)
      fr2(ir)=aimag(z1)
    end do
    call fderiv(-1,nri,r,fr1,gr)
    t1=gr(nri)
    call fderiv(-1,nri,r,fr2,gr)
    zfmtinp=(fourpi/dble(lmmaxinr))*cmplx(t1,gr(nri),8)
  else
    zfmtinp=0.d0
  end if
  nr0=nr-nri
  if (nr0.eq.0) return
  ir0=nri+1
  do ir=ir0,nr
    z1=zdotc(lmmaxvr,zfmt1(:,ir),1,zfmt2(:,ir),1)*(r(ir)**2)
    fr1(ir)=dble(z1)
    fr2(ir)=aimag(z1)
  end do
  call fderiv(-1,nr0,r(ir0),fr1(ir0),gr(ir0))
  t1=gr(nr)
  call fderiv(-1,nr0,r(ir0),fr2(ir0),gr(ir0))
  zfmtinp=zfmtinp+(fourpi/dble(lmmaxvr))*cmplx(t1,gr(nr),8)
end if
return
end function
!EOC


