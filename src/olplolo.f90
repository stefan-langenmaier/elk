
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olplolo(ias,ngp,o)
use modmain
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp
complex(8), intent(inout) :: o(*)
! local variables
integer ld,is,ilo,jlo
integer l,m,lm,i,j,k
ld=ngp+nlotot
is=idxis(ias)
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do jlo=1,nlorb(is)
    if (lorbl(jlo,is).eq.l) then
      do m=-l,l
        lm=idxlm(l,m)
        i=ngp+idxlo(lm,ilo,ias)
        j=ngp+idxlo(lm,jlo,ias)
        if (i.le.j) then
          k=i+(j-1)*ld
          o(k)=o(k)+ololo(ilo,jlo,ias)
        end if
      end do
    end if
  end do
end do
return
end subroutine

