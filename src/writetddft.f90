
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetddft
use modmain
use modtddft
use modmpi
implicit none
! local variables
character(256) fext
if (.not.mp_mpi) goto 10
! delete all files at first time-step
if (itimes.le.1) then
  open(50,file='MOMENT_TD.OUT')
  close(50,status='DELETE')
  open(50,file='MOMENTM_TD.OUT')
  close(50,status='DELETE')
end if
! write non-optional quantities
open(50,file='MOMENT_TD.OUT',action='WRITE',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),momtot(1:ndmag)
close(50)
open(50,file='MOMENTM_TD.OUT',action='WRITE',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),momtotm
close(50)
!**********

! write optional quantities
if (ntwrite.le.0) goto 10
if (mod(itimes-1,ntwrite).ne.0) goto 10
! file extension
write(fext,'("_TS",I8.8,".OUT")') itimes
! charge density in 1D
if (tdrho1d) then
  open(50,file='RHO1D'//trim(fext),action='WRITE',form='FORMATTED')
  open(51,file='RHOLINES'//trim(fext),action='WRITE',form='FORMATTED')
  call plot1d(50,51,1,rhomt,rhoir)
  close(50)
  close(51)
end if
!*******

10 continue
! synchronise MPI processes
call mpi_barrier(mpi_comm_kpt,ierror)
return
end subroutine

