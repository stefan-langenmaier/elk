
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdelete
use modmain
use modphonon
use modmpi
implicit none
character(256) fext
! only the master process should delete the files
if (mp_mpi) then
! construct the phonon file extension
  call phfext(iqph,isph,iaph,ipph,fext)
! delete the eigenvalue/vector files
  open(70,file='DEVALFV'//trim(fext))
  close(70,status='DELETE')
  open(70,file='DEVECFV'//trim(fext))
  close(70,status='DELETE')
  open(70,file='DEVECSV'//trim(fext))
  close(70,status='DELETE')
end if
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
return
end subroutine

