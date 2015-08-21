
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaa(ias,ngp,apwalm,o)
use modmain
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(inout) :: o(*)
! local variables
integer ld,is,l,m,lm,io
ld=ngp+nlotot
is=idxis(ias)
do l=0,lmaxmat
  do m=-l,l
    lm=idxlm(l,m)
    do io=1,apword(l,is)
      call zher2a(ngp,1.d0,apwalm(:,io,lm,ias),ld,o)
    end do
  end do
end do
return
end subroutine

