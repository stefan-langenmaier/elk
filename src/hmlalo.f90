
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlalo(tapp,is,ia,ngp,apwalm,v,h)
use modmain
implicit none
! arguments
logical, intent(in) :: tapp
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: v(*)
complex(8), intent(inout) :: h(*)
! local variables
integer ias,ld,io,ilo
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3
integer ist,i,j,k,ki,kj
complex(8) zsum
! allocatable arrays
complex(8), allocatable :: hi(:)
if (tapp) allocate(hi(ngp))
ld=ngp+nlotot
ias=idxas(ia,is)
do ilo=1,nlorb(is)
  l1=lorbl(ilo,is)
  do m1=-l1,l1
    lm1=idxlm(l1,m1)
    i=ngp+idxlo(lm1,ilo,ias)
    if (tapp) hi(:)=0.d0
    do l3=0,lmaxmat
      do m3=-l3,l3
        lm3=idxlm(l3,m3)
        do io=1,apword(l3,is)
          zsum=0.d0
          do l2=0,lmaxvr
            if (mod(l1+l2+l3,2).eq.0) then
              do m2=-l2,l2
                lm2=idxlm(l2,m2)
                zsum=zsum+gntyry(lm1,lm2,lm3)*hloa(ilo,io,l3,lm2,ias)
              end do
            end if
          end do
! note that what is actually computed is the Hermitian conjugate of <lo|H|APW>
          if (abs(dble(zsum))+abs(aimag(zsum)).gt.1.d-14) then
            if (tapp) then
! calculate the row of matrix H for application to v
              hi(1:ngp)=hi(1:ngp)+zsum*apwalm(1:ngp,io,lm3,ias)
            else
! calculate the matrix elements
              if (tpmat) then
                k=((i-1)*i)/2
              else
                k=(i-1)*ld
              end if
              do j=1,ngp
                k=k+1
                h(k)=h(k)+conjg(zsum*apwalm(j,io,lm3,ias))
              end do
            end if
          end if
        end do
      end do
    end do
! apply row of H to v
    if (tapp) then
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(j,k,ki,kj)
!$OMP DO
      do ist=1,nstfv
        k=(ist-1)*nmatmax
        ki=k+i
        do j=1,ngp
          kj=k+j
          h(ki)=h(ki)+hi(j)*v(kj)
          h(kj)=h(kj)+conjg(hi(j))*v(ki)
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
    end if
  end do
end do
if (tapp) deallocate(hi)
return
end subroutine

