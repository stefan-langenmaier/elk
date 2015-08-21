
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfmtadd(nr,nri,za,zfmt1,zfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: za
complex(8), intent(in) :: zfmt1(lmmaxvr,nr)
complex(8), intent(inout) :: zfmt2(lmmaxvr,nr)
! local variables
integer nr0,ir0,ir
! add on inner part of muffin-tin
do ir=1,nri
  call zaxpy(lmmaxinr,za,zfmt1(:,ir),1,zfmt2(:,ir),1)
end do
! add on outer part of muffin-tin
nr0=nr-nri
if (nr0.eq.0) return
ir0=nri+1
call zaxpy(lmmaxvr*nr0,za,zfmt1(:,ir0),1,zfmt2(:,ir0),1)
return
end subroutine

