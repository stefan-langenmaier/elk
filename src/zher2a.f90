
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zher2a(n,alpha,x,ld,a)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: alpha
complex(8), intent(in) :: x(n)
integer, intent(in) :: ld
complex(8), intent(inout) :: a(*)
! local variables
integer j,k
! numbers less than eps are considered to be zero
real(8), parameter :: eps=1.d-10
real(8) a1,b1
complex(8) zt1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,zt1,a1,b1)
!$OMP DO
do j=1,n
  k=(j-1)*ld
  zt1=alpha*x(j)
  if (abs(dble(zt1)).gt.eps) then
    if (abs(aimag(zt1)).gt.eps) then
! complex prefactor
      a(k+1:k+j-1)=a(k+1:k+j-1)+zt1*conjg(x(1:j-1))
      a(k+j)=dble(a(k+j))+dble(zt1*conjg(x(j)))
    else
! real prefactor
      a1=dble(zt1)
      a(k+1:k+j-1)=a(k+1:k+j-1)+a1*conjg(x(1:j-1))
      a(k+j)=dble(a(k+j))+a1*dble(x(j))
    end if
  else
    if (abs(aimag(zt1)).gt.eps) then
! imaginary prefactor
      b1=aimag(zt1)
      a(k+1:k+j-1)=a(k+1:k+j-1)+b1*cmplx(aimag(x(1:j-1)),dble(x(1:j-1)),8)
      a(k+j)=dble(a(k+j))+b1*aimag(x(j))
    end if
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

