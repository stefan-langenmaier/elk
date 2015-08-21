
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mae
use modmain
use modmpi
use modstore
implicit none
! local variables
integer itp,im(2)
real(8), parameter :: bfm=2.d0
real(8) em(2),de
! automatic arrays
real(8) tp(2,npmae)
! store original parameters
spinpol0=spinpol
spinorb0=spinorb
bfieldc00(:)=bfieldc0(:)
reducebf0=reducebf
fixspin0=fixspin
momfix0(:)=momfix(:)
! generate a spherical cover
call sphcover(npmae,tp)
! initialise from atomic densities for all ground-state runs
trdstate=.false.
! spin-orbit coupling must be enabled
spinorb=.true.
! reduce the external magnetic field after each s.c. loop
reducebf=0.75d0
! global fixed spin moment direction
fixspin=-1
! open MAE_INFO.OUT
if (mp_mpi) open(71,file='MAE_INFO.OUT',action='WRITE',form='FORMATTED')
im(:)=1
em(1)=1.d8
em(2)=-1.d8
! loop over points on sphere
do itp=1,npmae
  if (mp_mpi) then
    write(*,'("Info(findmae): fixed spin moment direction ",I6," of ",I6)') &
     itp,npmae
  end if
! set fixed spin moment direction
  momfix(1)=sin(tp(1,itp))*cos(tp(2,itp))
  momfix(2)=sin(tp(1,itp))*sin(tp(2,itp))
  momfix(3)=cos(tp(1,itp))
! large magnetic field in the opposite direction as fixed moment
  bfieldc0(:)=-bfm*momfix(:)
! run the ground-state calculation
  call gndstate
  if (mp_mpi) then
    write(71,*)
    write(71,'("Fixed spin moment direction point ",I6," of ",I6)') itp,npmae
    write(71,'("Spherical coordinates of direction : ",2G18.10)') tp(:,itp)
    write(71,'("Direction vector : ",3G18.10)') momfix(:)
    write(71,'("Calculated total moment magnitude : ",G18.10)') momtotm
    write(71,'("Total energy : ",G22.12)') engytot
    call flushifc(71)
  end if
! check for minimum and maximum total energy
  if (engytot.lt.em(1)) then
    em(1)=engytot
    im(1)=itp
  end if
  if (engytot.gt.em(2)) then
    em(2)=engytot
    im(2)=itp
  end if
end do
! magnetic anisotropy energy
de=em(2)-em(1)
if (mp_mpi) then
  write(71,*)
  write(71,'("Minimum energy point : ",I6)') im(1)
  write(71,'("Maximum energy point : ",I6)') im(2)
  write(71,*)
  write(71,'("Estimated magnetic anisotropy energy (MAE) : ",G18.10)') de
  write(71,*)
  write(71,'("MAE per unit volume : ",G18.10)') de/omega
  close(71)
  open(50,file='MAE.OUT',action='WRITE',form='FORMATTED')
  write(50,'(G18.10)') de
  close(50)
  open(50,file='MAEPUV.OUT',action='WRITE',form='FORMATTED')
  write(50,'(G18.10)') de/omega
  close(50)
  write(*,*)
  write(*,'("Info(mae):")')
  write(*,'(" Estimated magnetic anisotropy energy written to MAE.OUT")')
  write(*,'(" MAE per unit volume written to MAEPUV.OUT")')
  write(*,*)
  write(*,'(" Number of fixed spin moment directions used : ",I6)') npmae
  write(*,*)
  write(*,'(" Additional information written to MAE_INFO.OUT")')
end if
! restore original input parameters
spinpol=spinpol0
spinorb=spinorb0
bfieldc0(:)=bfieldc00(:)
reducebf=reducebf0
fixspin=fixspin0
momfix(:)=momfix0(:)
return
end subroutine

