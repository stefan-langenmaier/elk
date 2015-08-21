
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dolpaa(ias,ngp,ngpq,apwalm,apwalmq,dapwalm,dapwalmq,ld,od)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp,ngpq
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw)
complex(8), intent(in) :: dapwalmq(ngkmax,apwordmax,lmmaxapw)
integer, intent(in) :: ld
complex(8), intent(inout) :: od(ld,*)
! local variables
integer is,l,m,lm,io
if (ias.ne.iasph) return
is=idxis(ias)
lm=0
do l=0,lmaxmat
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      call zgerci(ngpq,ngp,zone,apwalmq(:,io,lm,ias),dapwalm(:,io,lm),ld,od)
      call zgerci(ngpq,ngp,zone,dapwalmq(:,io,lm),apwalm(:,io,lm,ias),ld,od)
    end do
  end do
end do
return
end subroutine

