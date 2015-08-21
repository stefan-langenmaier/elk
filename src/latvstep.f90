
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine latvstep
use modmain
use modmpi
use modstore
implicit none
integer i
real(8) t1
do i=1,3
! compute the dot-product between the current and previous stress vector
  t1=dot_product(stress(:,i),stressp(:,i))
! if stress vector is in the same direction then increase step size parameter
  if (t1.gt.0.d0) then
    taulatv(i)=taulatv(i)+tau0latv
  else
    taulatv(i)=tau0latv
  end if
  avec(:,i)=avec(:,i)-taulatv(i)*(stress(:,i)+stressp(:,i))
end do
! scale the vectors to conserve volume if required
if (latvopt.eq.2) then
  call reciplat
  t1=(omega0/omega)**(1.d0/3.d0)
  avec(:,:)=t1*avec(:,:)
end if
! each MPI process should have identical lattice vectors
call mpi_bcast(avec,9,mpi_double_precision,0,mpi_comm_kpt,ierror)
return
end subroutine

