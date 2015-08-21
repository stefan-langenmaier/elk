
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine geomopt
use modmain
use modmpi
implicit none
! local variables
integer istp
! initialise global variables
call init0
! muffin-tin radii must be fixed during optimisation
autormt=.false.
! atomic forces are required
tforce=.true.
if (task.eq.3) then
  trdstate=.true.
else
  trdstate=.false.
end if
! initial step sizes
if (allocated(tauatm)) deallocate(tauatm)
allocate(tauatm(natmtot))
tauatm(:)=tau0atm
! initialise the previous total force on each atom
if (allocated(forcetotp)) deallocate(forcetotp)
allocate(forcetotp(3,natmtot))
forcetotp(:,:)=0.d0
! open TOTENERGY.OUT
open(71,file='TOTENERGY_OPT.OUT',action='WRITE',form='FORMATTED')
! open FORCEMAX.OUT
open(72,file='FORCEMAX.OUT',action='WRITE',form='FORMATTED')
! open GEOMETRY_OPT.OUT
open(73,file='GEOMETRY_OPT.OUT',action='WRITE',form='FORMATTED')
! open IADIST_OPT.OUT
open(74,file='IADIST_OPT.OUT',action='WRITE',form='FORMATTED')
! open FORCES_OPT.OUT
open(75,file='FORCES_OPT.OUT',action='WRITE',form='FORMATTED')
do istp=1,maxgeostp
  if (mp_mpi) write(*,'("Info(geomopt): geometry optimisation step : ",I6)') &
   istp
! ground-state and forces calculation
  call gndstate
! check that calculation converged
  if (iscl.ge.maxscl) then
    write(*,*)
    write(*,'("Warning(geomopt): ground-state calculation failed to converge")')
  end if
! subsequent calculations will read in the potential from STATE.OUT
  trdstate=.true.
  if (mp_mpi) then
! write the converged total energy to TOTENERGY_OPT.OUT
    write(71,'(G22.12)') engytot
    call flushifc(71)
! write maximum force magnitude to FORCEMAX.OUT
    write(72,'(G18.10)') forcemax
    call flushifc(72)
  end if
! update the atomic positions
  call geomstep
! write optimised atomic positions, interatomic distances and forces to file
  if (mp_mpi) then
    write(73,*); write(73,*)
    write(73,'("! Geometry optimisation step : ",I6)') istp
    call writegeom(73)
    call flushifc(73)
    write(74,*); write(74,*)
    write(74,'("Geometry optimisation step : ",I6)') istp
    call writeiad(74)
    call flushifc(74)
    write(75,*); write(75,*)
    write(75,'("Geometry optimisation step : ",I6)') istp
    call writeforces(75)
    write(75,*)
    write(75,'("Maximum force magnitude over all atoms (target) : ",G18.10,&
     &" (",G18.10,")")') forcemax,epsforce
    call flushifc(75)
  end if
! check force convergence
  if (forcemax.le.epsforce) then
    write(75,*)
    write(75,'("Force convergence target achieved")')
    goto 10
  end if
end do
write(*,*)
write(*,'("Warning(geomopt): geometry optimisation failed to converge in ",I6,&
 &" steps")') maxgeostp
10 continue
close(71); close(72); close(73); close(74); close(75)
return
end subroutine

