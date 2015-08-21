
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rbsht(nr,nri,rfmt1,rfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: rfmt1(lmmaxvr,nr)
real(8), intent(out) :: rfmt2(lmmaxvr,nr)
! local variables
integer nro,ir
! transform the inner part of the muffin-tin
call dgemm('N','N',lmmaxvr,nri,lmmaxinr,1.d0,rbshtvr,lmmaxvr,rfmt1,lmmaxvr, &
 0.d0,rfmt2,lmmaxvr)
! transform the outer part of the muffin-tin
if (nr.eq.nri) return
nro=nr-nri
ir=nri+1
call dgemm('N','N',lmmaxvr,nro,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rfmt1(:,ir), &
 lmmaxvr,0.d0,rfmt2(:,ir),lmmaxvr)
return
end subroutine

