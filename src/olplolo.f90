
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olplolo(tapp,is,ia,ngp,v,o)
use modmain
implicit none
! arguments
logical, intent(in) :: tapp
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: v(*)
complex(8), intent(inout) :: o(*)
! local variables
integer ld,ias,ilo1,ilo2,l,m,lm
integer ist,i,j,k,ki,kj
ld=ngp+nlotot
ias=idxas(ia,is)
do ilo1=1,nlorb(is)
  l=lorbl(ilo1,is)
  do ilo2=1,nlorb(is)
    if (lorbl(ilo2,is).eq.l) then
      do m=-l,l
        lm=idxlm(l,m)
        i=ngp+idxlo(lm,ilo1,ias)
        j=ngp+idxlo(lm,ilo2,ias)
        if (i.le.j) then
          if (tapp) then
! apply the overlap operator to v
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,ki,kj)
!$OMP DO
            do ist=1,nstfv
              k=(ist-1)*nmatmax
              ki=k+i
              kj=k+j
              o(ki)=o(ki)+ololo(ilo1,ilo2,ias)*v(kj)
              if (i.ne.j) o(kj)=o(kj)+ololo(ilo1,ilo2,ias)*v(ki)
            end do
!$OMP END DO
!$OMP END PARALLEL
          else
! calculate the matrix elements
            if (tpmat) then
              k=i+((j-1)*j)/2
            else
              k=i+(j-1)*ld
            end if
            o(k)=o(k)+ololo(ilo1,ilo2,ias)
          end if
        end if
      end do
    end if
  end do
end do
return
end subroutine

