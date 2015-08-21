
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfmtinp
! !INTERFACE:
complex(8) function zfmtinp(nr,nri,r,zfmt1,zfmt2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of points on the inner part of the muffin-tin (in,integer)
!   r     : radial mesh (in,real(nr))
!   zfmt1 : first complex muffin-tin function in spherical harmonics
!           (in,complex(lmmaxvr,nr))
!   zfmt2 : second complex muffin-tin function (in,complex(lmmaxvr,nr))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions in the muffin-tin. In
!   other words, given two complex functions of the form
!   $$ f({\bf r})=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}f_{lm}(r)Y_{lm}
!    (\hat{\bf r}), $$
!   the function returns
!   $$ I=\sum_{l=0}^{l_{\rm max}}\sum_{m=-l}^{l}\int f_{lm}^{1*}(r)
!    f_{lm}^2(r)r^2\,dr\;. $$
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!   Modified, September 2013 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr)
complex(8), intent(in) :: zfmt1(lmmaxvr,nr),zfmt2(lmmaxvr,nr)
! local variables
integer ir
real(8) t1
complex(8) z1
! automatic arrays
real(8) fr1(nr),fr2(nr),gr(nr)
! external functions
complex(8) zdotc
external zdotc
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
return
end function
!EOC


