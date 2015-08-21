
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genstrain
use modmain
implicit none
! local variables
integer i,j,k,l
real(8) e(3,3,12),a(3,3),t1
! external functions
real(8) ddot,dnrm2
external ddot,dnrm2
! store rotations around x, y and z
e(:,:,:)=0.d0
t1=1.d0/sqrt(2.d0)
e(2,3,1)=t1
e(3,2,1)=-t1
e(1,3,2)=t1
e(3,1,2)=-t1
e(1,2,3)=t1
e(2,1,3)=-t1
k=3
nstrain=0
do i=1,3
  do j=1,3
! set strain tensor to delta_ij
    a(:,:)=0.d0
    a(i,j)=1.d0
! symmetrise strain tensor
    call symmat(a)
! orthogonalise strain tensor to previous tensors
    do l=1,k
      t1=-ddot(9,e(:,:,l),1,a,1)
      call daxpy(9,t1,e(:,:,l),1,a,1)
    end do
! compute the norm
    t1=dnrm2(9,a,1)
    if (t1.lt.epslat) cycle
! normalise tensor
    k=k+1
    e(:,:,k)=a(:,:)/t1
! store the strain tensor in global array
    nstrain=k-3
    strain(:,:,nstrain)=e(:,:,k)
  end do
end do
! zero small components
do k=1,nstrain
  do i=1,3
    do j=1,3
      if (abs(strain(i,j,k)).lt.epslat) strain(i,j,k)=0.d0
    end do
  end do
end do
return
end subroutine

