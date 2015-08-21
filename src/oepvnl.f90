
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepvnl(vnlcv,vnlvv)
use modmain
use modmpi
implicit none
! arguments
complex(8), intent(out) :: vnlcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(out) :: vnlvv(nstsv,nstsv,nkpt)
! local variables
integer ik,ncv,nvv,lp
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL
  write(*,'("Info(oepvnl): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
  call oepvnlk(ik,vnlcv(:,:,:,ik),vnlvv(:,:,ik))
end do
!$OMP END DO
!$OMP END PARALLEL
! broadcast matrix elements to all other processes
ncv=ncrmax*natmtot*nstsv
nvv=nstsv*nstsv
do ik=1,nkpt
  lp=mod(ik-1,np_mpi)
  call mpi_bcast(vnlcv(:,:,:,ik),ncv,mpi_double_complex,lp,mpi_comm_world, &
   ierror)
  call mpi_bcast(vnlvv(:,:,ik),nvv,mpi_double_complex,lp,mpi_comm_world,ierror)
end do
return
end subroutine

