
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zvmul1(n,a,b,c)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: a(n),b(n)
complex(8), intent(out) :: c(n)
c(:)=conjg(a(:))*b(:)
return
end subroutine

