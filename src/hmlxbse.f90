
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlxbse
use modmain
use modmpi
implicit none
! local variables
integer ik2
! allocatable arrays
real(8), allocatable :: jlgr(:,:,:)
! allocate local arrays
allocate(jlgr(0:lnpsd+1,ngvec,nspecies))
! compute the required spherical Bessel functions
call genjlgpr(lnpsd+1,gc,jlgr)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik2=1,nkptnr
! distribute among MPI processes
  if (mod(ik2-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL
  write(*,'("Info(hmlxbse): ",I6," of ",I6," k-points")') ik2,nkptnr
!$OMP END CRITICAL
  call hmlxbsek(ik2,jlgr)
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(jlgr)
return
end subroutine

