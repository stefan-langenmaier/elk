
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfmtmul2(nr,nri,x1,x2,y1,y2,z)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: x1(lmmaxvr,nr),x2(lmmaxvr,nr)
complex(8), intent(in) :: y1(lmmaxvr,nr),y2(lmmaxvr,nr)
complex(8), intent(out) :: z(lmmaxvr,nr)
! local variables
integer ir
do ir=1,nri
  z(1:lmmaxinr,ir)=conjg(x1(1:lmmaxinr,ir))*y1(1:lmmaxinr,ir) &
                  +conjg(x2(1:lmmaxinr,ir))*y2(1:lmmaxinr,ir)
end do
do ir=nri+1,nr
  z(:,ir)=conjg(x1(:,ir))*y1(:,ir)+conjg(x2(:,ir))*y2(:,ir)
end do
return
end subroutine

