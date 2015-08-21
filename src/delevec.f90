
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: delevec
! !INTERFACE:
subroutine delevec
! !USES:
use modmain
use modmpi
! !DESCRIPTION:
!   Deletes the first- and second-variational eigenvector files {\tt EVECFV.OUT}
!   and {\tt EVECSV.OUT}.
!
! !REVISION HISTORY:
!   Created May 2007 (JKD)
!EOP
!BOC
implicit none
! only the master process should delete the files
if (mp_mpi) then
! delete the first-variational eigenvector file
  open(70,file=trim(scrpath)//'EVECFV'//trim(filext))
  close(70,status='DELETE')
! delete the second-variational eigenvector file
  open(70,file=trim(scrpath)//'EVECSV'//trim(filext))
  close(70,status='DELETE')
end if
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
return
end subroutine
!EOC

