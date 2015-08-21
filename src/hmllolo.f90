
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmllolo(tapp,is,ia,ngp,v,h)
use modmain
implicit none
! arguments
logical, intent(in) :: tapp
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: v(*)
complex(8), intent(inout) :: h(*)
! local variables
integer ld,ias,ilo1,ilo2
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3
integer ist,i,j,k,ki,kj
complex(8) zsum
ld=ngp+nlotot
ias=idxas(ia,is)
do ilo1=1,nlorb(is)
  l1=lorbl(ilo1,is)
  do m1=-l1,l1
    lm1=idxlm(l1,m1)
    i=ngp+idxlo(lm1,ilo1,ias)
    do ilo2=1,nlorb(is)
      l3=lorbl(ilo2,is)
      do m3=-l3,l3
        lm3=idxlm(l3,m3)
        j=ngp+idxlo(lm3,ilo2,ias)
        if (i.le.j) then
          zsum=0.d0
          do l2=0,lmaxvr
            if (mod(l1+l2+l3,2).eq.0) then
              do m2=-l2,l2
                lm2=idxlm(l2,m2)
                zsum=zsum+gntyry(lm1,lm2,lm3)*hlolo(ilo1,ilo2,lm2,ias)
              end do
            end if
          end do
          if (tapp) then
! apply the Hamiltonian operator to v
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(k,ki,kj)
!$OMP DO
            do ist=1,nstfv
              k=(ist-1)*nmatmax
              ki=k+i
              kj=k+j
              h(ki)=h(ki)+zsum*v(kj)
              if (i.ne.j) h(kj)=h(kj)+conjg(zsum)*v(ki)
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
            h(k)=h(k)+zsum
          end if
        end if
      end do
    end do
  end do
end do
return
end subroutine

