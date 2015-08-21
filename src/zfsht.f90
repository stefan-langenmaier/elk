
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfsht(nr,nri,zfmt1,zfmt2)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
complex(8), intent(in) :: zfmt1(lmmaxvr,nr)
complex(8), intent(out) :: zfmt2(lmmaxvr,nr)
! local variables
integer nr0,ir0
! transform the inner part of the muffin-tin
call zgemm('N','N',lmmaxinr,nri,lmmaxinr,zone,zfshtinr,lmmaxinr,zfmt1,lmmaxvr, &
 zzero,zfmt2,lmmaxvr)
! transform the outer part of the muffin-tin
nr0=nr-nri
if (nr0.eq.0) return
ir0=nri+1
call zgemm('N','N',lmmaxvr,nr0,lmmaxvr,zone,zfshtvr,lmmaxvr,zfmt1(:,ir0), &
 lmmaxvr,zzero,zfmt2(:,ir0),lmmaxvr)
return
end subroutine

