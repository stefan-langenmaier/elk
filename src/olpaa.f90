
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpaa(tapp,is,ia,ngp,apwalm,v,o)
use modmain
implicit none
! arguments
logical, intent(in) :: tapp
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: v(*)
complex(8), intent(inout) :: o(*)
! local variables
integer ld,ias,l,m,lm,io
ld=ngp+nlotot
ias=idxas(ia,is)
do l=0,lmaxmat
  do m=-l,l
    lm=idxlm(l,m)
    do io=1,apword(l,is)
      if (tapp) then
! apply the overlap to a set of vectors
        call zmatinpv(ngp,zhalf,apwalm(:,io,lm,ias),apwalm(:,io,lm,ias),nstfv, &
         nmatmax,v,o)
      else
! compute the matrix explicitly
        call zmatinp(tpmat,ngp,zhalf,apwalm(:,io,lm,ias),apwalm(:,io,lm,ias), &
         ld,o)
      end if
    end do
  end do
end do
return
end subroutine

