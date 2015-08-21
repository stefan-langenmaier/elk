
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zvmul2(n,a1,b1,a2,b2,c)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: a1(n),b1(n),a2(n),b2(n)
complex(8), intent(out) :: c(n)
c(:)=conjg(a1(:))*b1(:)+conjg(a2(:))*b2(:)
return
end subroutine

