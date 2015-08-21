
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symmat(ac)
use modmain
implicit none
! arguments
real(8), intent(inout) :: ac(3,3)
! local variables
integer isym,lsp
real(8) as(3,3),b(3,3),c(3,3),t1
as(:,:)=0.d0
do isym=1,nsymcrys
  lsp=lsplsymc(isym)
  call r3mm(symlatc(:,:,lsp),ac,b)
  call r3mmt(b,symlatc(:,:,lsp),c)
  as(:,:)=as(:,:)+c(:,:)
end do
t1=1.d0/dble(nsymcrys)
ac(:,:)=t1*as(:,:)
return
end subroutine

