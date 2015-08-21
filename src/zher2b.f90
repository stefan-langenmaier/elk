
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: zher2b
! !INTERFACE:
subroutine zher2b(n,alpha,x,y,ld,a)
! !INPUT/OUTPUT PARAMETERS:
!   n     : length of vectors (in,integer)
!   alpha : complex constant (in,complex)
!   x     : first input vector (in,complex(n))
!   y     : second input vector (in,complex(n))
!   ld    : leading dimension of A, not used if A is packed (in,integer)
!   a     : output matrix A (inout,complex(*))
! !DESCRIPTION:
!   Performs the rank-2 operation
!   $$ A_{ij}\rightarrow\alpha{\bf x}_i^*{\bf y}_j+\alpha^*{\bf y}_i^*{\bf x}_j
!    +A_{ij}. $$
!   This is similar to the {\tt BLAS} routine {\tt zher2}, except that here a
!   matrix of inner products is formed instead of an outer product of vectors.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: alpha
complex(8), intent(in) :: x(n)
complex(8), intent(in) :: y(n)
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
  zt1=alpha*y(j)
  if (abs(dble(zt1)).gt.eps) then
    if (abs(aimag(zt1)).gt.eps) then
! complex prefactor
      a(k+1:k+j-1)=a(k+1:k+j-1)+zt1*conjg(x(1:j-1))
      a(k+j)=dble(a(k+j))+2.d0*dble(zt1*conjg(x(j)))
    else
! real prefactor
      a1=dble(zt1)
      a(k+1:k+j-1)=a(k+1:k+j-1)+a1*conjg(x(1:j-1))
      a(k+j)=dble(a(k+j))+2.d0*a1*dble(x(j))
    end if
  else
    if (abs(aimag(zt1)).gt.eps) then
! imaginary prefactor
      b1=aimag(zt1)
      a(k+1:k+j-1)=a(k+1:k+j-1)+b1*cmplx(aimag(x(1:j-1)),dble(x(1:j-1)),8)
      a(k+j)=dble(a(k+j))+2.d0*b1*aimag(x(j))
    end if
  end if
  zt1=conjg(alpha)*x(j)
  if (abs(dble(zt1)).gt.eps) then
    if (abs(aimag(zt1)).gt.eps) then
! complex prefactor
      a(k+1:k+j-1)=a(k+1:k+j-1)+zt1*conjg(y(1:j-1))
    else
! real prefactor
      a1=dble(zt1)
      a(k+1:k+j-1)=a(k+1:k+j-1)+a1*conjg(y(1:j-1))
    end if
  else
    if (abs(aimag(zt1)).gt.eps) then
! imaginary prefactor
      b1=aimag(zt1)
      a(k+1:k+j-1)=a(k+1:k+j-1)+b1*cmplx(aimag(y(1:j-1)),dble(y(1:j-1)),8)
    end if
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC
