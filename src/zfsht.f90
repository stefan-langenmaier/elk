
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
integer nro,ir
! transform the inner part of the muffin-tin
call zgemm('N','N',lmmaxinr,nri,lmmaxvr,zone,zfshtvr,lmmaxvr,zfmt1,lmmaxvr, &
 zzero,zfmt2,lmmaxvr)
zfmt2(lmmaxinr+1:lmmaxvr,1:nri)=0.d0
! transform the outer part of the muffin-tin
if (nr.eq.nri) return
nro=nr-nri
ir=nri+1
call zgemm('N','N',lmmaxvr,nro,lmmaxvr,zone,zfshtvr,lmmaxvr,zfmt1(:,ir), &
 lmmaxvr,zzero,zfmt2(:,ir),lmmaxvr)
return
end subroutine

