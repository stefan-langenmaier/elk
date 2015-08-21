
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writevnlijji
! !INTERFACE:
subroutine writevnlijji
! !USES:
use modmain
use modmpi
! !DESCRIPTION:
!   Generates non-local Coulomb matrix elements of the type $(i-jj-i)$ and
!   outputs them to the file {\tt VNLIJJI.OUT}.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! allocatable arrays
real(8), allocatable :: vnlijji(:,:,:)
integer recl,ik
! determine record length for vnlijji and open file
allocate(vnlijji(nstsv,nstsv,nkpt))
inquire(iolength=recl) vnlijji
deallocate(vnlijji)
open(100,file='VNLIJJI.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vnlijji)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(vnlijji(nstsv,nstsv,nkpt))
!$OMP CRITICAL
  write(*,'("Info(writevnlijji): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
! calculate non-local matrix elements of the type (i-jj-i)
  call genvnlijji(ik,vnlijji)
!$OMP CRITICAL
  write(100,rec=ik) vnlijji
!$OMP END CRITICAL
  deallocate(vnlijji)
end do
!$OMP END DO
!$OMP END PARALLEL
close(100)
return
end subroutine
!EOC
