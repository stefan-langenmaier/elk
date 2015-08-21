
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writestress
use modmain
use modmpi
implicit none
! local variables
integer i,j
! initialise universal variables
call init0
! start from the atomic densities
trdstate=.false.
! generate the stress matrix
call genstress
! write the stress matrix to file
if (mp_mpi) then
  open(50,file='STRESS.OUT',action='WRITE',form='FORMATTED')
  do j=1,3
    write(50,*)
    write(50,'("Derivative of total energy w.r.t. lattice vector ",I1," :")') j
    write(50,'(3G18.10)') (stress(i,j),i=1,3)
  end do
  close(50)
  write(*,*)
  write(*,'("Info(writestress):")')
  write(*,'(" Stress matrix written to STRESS.OUT")')
end if
return
end subroutine

