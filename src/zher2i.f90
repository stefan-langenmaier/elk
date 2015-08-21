
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zher2i(n,alpha,x,y,ld,a)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: alpha
complex(8), intent(in) :: x(n),y(n)
integer, intent(in) :: ld
complex(8), intent(inout) :: a(*)
! local variables
integer j,k
! numbers less than eps are considered to be zero
real(8), parameter :: eps=1.d-10
real(8) a1,b1
complex(8) z1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,z1,a1,b1)
!$OMP DO
do j=1,n
  k=(j-1)*ld
  z1=alpha*y(j)
  if (abs(dble(z1)).gt.eps) then
    if (abs(aimag(z1)).gt.eps) then
! complex prefactor
      a(k+1:k+j-1)=a(k+1:k+j-1)+z1*conjg(x(1:j-1))
      a(k+j)=dble(a(k+j))+dble(z1*conjg(x(j)))
    else
! real prefactor
      a1=dble(z1)
      a(k+1:k+j-1)=a(k+1:k+j-1)+a1*conjg(x(1:j-1))
      a(k+j)=dble(a(k+j))+a1*dble(x(j))
    end if
  else
    if (abs(aimag(z1)).gt.eps) then
! imaginary prefactor
      b1=aimag(z1)
      a(k+1:k+j-1)=a(k+1:k+j-1)+b1*cmplx(aimag(x(1:j-1)),dble(x(1:j-1)),8)
      a(k+j)=dble(a(k+j))+b1*aimag(x(j))
    end if
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

