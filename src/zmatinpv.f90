
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zmatinpv(n,alpha,x,y,nv,ld,v,av)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: alpha
complex(8), intent(in) :: x(n)
complex(8), intent(in) :: y(n)
integer, intent(in) :: nv
integer, intent(in) :: ld
complex(8), intent(in) :: v(ld,nv)
complex(8), intent(inout) :: av(ld,nv)
! local variables
integer j
! numbers less than eps are considered to be zero
real(8), parameter :: eps=1.d-10
real(8) a1,b1
complex(8) zt1
! external functions
complex(8) zdotu
external zdotu
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zt1,a1,b1)
!$OMP DO
do j=1,nv
  zt1=zdotu(n,x,1,v(:,j),1)
  zt1=alpha*conjg(zt1)
  if (abs(dble(zt1)).gt.eps) then
    if (abs(aimag(zt1)).gt.eps) then
! complex prefactor
      av(1:n,j)=av(1:n,j)+conjg(zt1*y(1:n))
    else
! real prefactor
      a1=dble(zt1)
      av(1:n,j)=av(1:n,j)+a1*conjg(y(1:n))
    end if
  else
! imaginary prefactor
    if (abs(aimag(zt1)).gt.eps) then
      b1=aimag(zt1)
      av(1:n,j)=av(1:n,j)-b1*cmplx(aimag(y(1:n)),dble(y(1:n)),8)
    end if
  end if
  zt1=zdotu(n,y,1,v(:,j),1)
  zt1=conjg(alpha*zt1)
  if (abs(dble(zt1)).gt.eps) then
    if (abs(aimag(zt1)).gt.eps) then
      av(1:n,j)=av(1:n,j)+conjg(zt1*x(1:n))
    else
      a1=dble(zt1)
      av(1:n,j)=av(1:n,j)+a1*conjg(x(1:n))
    end if
  else
    if (abs(aimag(zt1)).gt.eps) then
      b1=aimag(zt1)
      av(1:n,j)=av(1:n,j)-b1*cmplx(aimag(x(1:n)),dble(x(1:n)),8)
    end if
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
