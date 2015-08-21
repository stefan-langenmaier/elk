
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine olpalo(tapp,is,ia,ngp,apwalm,v,o)
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
integer ld,ias,ilo,io,l,m,lm
integer ist,i,j,k,ki,kj
! allocatable arrays
complex(8), allocatable :: oj(:)
if (tapp) allocate(oj(ngp))
ld=ngp+nlotot
ias=idxas(ia,is)
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do m=-l,l
    lm=idxlm(l,m)
    j=ngp+idxlo(lm,ilo,ias)
    if (tapp) then
! calculate the column of matrix elements for application to v
      do i=1,ngp
        oj(i)=0.d0
        do io=1,apword(l,is)
          oj(i)=oj(i)+conjg(apwalm(i,io,lm,ias))*oalo(io,ilo,ias)
        end do
      end do
! apply column of O to v
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,k,ki,kj)
!$OMP DO
      do ist=1,nstfv
        k=(ist-1)*nmatmax
        kj=k+j
        do i=1,ngp
          ki=k+i
          o(ki)=o(ki)+oj(i)*v(kj)
          o(kj)=o(kj)+conjg(oj(i))*v(ki)
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
    else
! calculate the matrix elements
      if (tpmat) then
        k=((j-1)*j)/2
      else
        k=(j-1)*ld
      end if
      do i=1,ngp
        k=k+1
        do io=1,apword(l,is)
          o(k)=o(k)+conjg(apwalm(i,io,lm,ias))*oalo(io,ilo,ias)
        end do
      end do
    end if
  end do
end do
if (tapp) deallocate(oj)
return
end subroutine

