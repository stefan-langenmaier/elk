
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfmtmul1(nr,nri,x,y,z)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: x(lmmaxvr,nr),y(lmmaxvr,nr)
complex(8), intent(out) :: z(lmmaxvr,nr)
! local variables
integer ir
do ir=1,nri
  z(1:lmmaxinr,ir)=conjg(x(1:lmmaxinr,ir))*y(1:lmmaxinr,ir)
end do
do ir=nri+1,nr
  z(:,ir)=conjg(x(:,ir))*y(:,ir)
end do
return
end subroutine

