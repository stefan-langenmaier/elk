
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genstress
use modmain
use modmpi
use modstore
implicit none
! local variables
integer i
real(8) a(3,3),et1,t1
! store original parameters
avec0(:,:)=avec(:,:)
tshift0=tshift
tshift=.false.
tforce0=tforce
tforce=.false.
! generate the strain tensor
call genstrain
! loop over strain tensors
do i=1,nstrain
  if (mp_mpi) then
    write(*,'("Info(genstress): strain tensor : ",I1)') i
  end if
! restore orginal lattice vectors and symmetries
  avec(:,:)=avec0(:,:)
  call symmetry
! displace lattice vectors
  call r3mm(strain(:,:,i),avec,a)
  avec(:,:)=avec(:,:)+deltast*a(:,:)
! run the ground-state calculation
  call gndstate
! subsequent calculations will read STATE.OUT
  trdstate=.true.
! store the total energy for the first displacement
  et1=engytot
! displace the lattice vector again
  avec(:,:)=avec(:,:)+deltast*a(:,:)
! run the ground-state calculation again
  call gndstate
! compute the stress tensor component
  stress(i)=(engytot-et1)/deltast
end do
! compute the maximum stress magnitude over all lattice vectors
stressmax=0.d0
do i=1,nstrain
  t1=abs(stress(i))
  if (t1.gt.stressmax) stressmax=t1
end do
! restore original parameters
avec(:,:)=avec0(:,:)
tshift=tshift0
tforce=tforce0
return
end subroutine

