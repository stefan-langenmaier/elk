
! Copyright (C) 2014 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function rzfmtinp(nr,nri,r,rfmt,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: r(nr)
real(8), intent(in) :: rfmt(lmmaxvr,nr)
complex(8), intent(in) :: zfmt(lmmaxvr,nr)
! local variables
integer ir,itp
real(8) t1
complex(8) z1
! automatic arrays
real(8) fr1(nr),fr2(nr),gr(nr)
t1=fourpi/dble(lmmaxinr)
do ir=1,nri
  z1=rfmt(1,ir)*zfmt(1,ir)
  do itp=2,lmmaxinr
    z1=z1+rfmt(itp,ir)*zfmt(itp,ir)
  end do
  z1=z1*(t1*r(ir)**2)
  fr1(ir)=dble(z1)
  fr2(ir)=aimag(z1)
end do
t1=fourpi/dble(lmmaxvr)
do ir=nri+1,nr
  z1=rfmt(1,ir)*zfmt(1,ir)
  do itp=2,lmmaxvr
    z1=z1+rfmt(itp,ir)*zfmt(itp,ir)
  end do
  z1=z1*(t1*r(ir)**2)
  fr1(ir)=dble(z1)
  fr2(ir)=aimag(z1)
end do
call fderiv(-1,nr,r,fr1,gr)
t1=gr(nr)
call fderiv(-1,nr,r,fr2,gr)
rzfmtinp=cmplx(t1,gr(nr),8)
return
end function

