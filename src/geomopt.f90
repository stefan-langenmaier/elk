
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine geomopt
use modmain
use modmpi
use modstore
implicit none
! local variables
integer istp,jstp,i,j
real(8) t1
! initialise global variables
call init0
! store orginal volume
omega0=omega
! atomic forces are required
tforce=.true.
if (task.eq.3) then
  trdstate=.true.
else
  trdstate=.false.
end if
! initial atomic step sizes
if (allocated(tauatp)) deallocate(tauatp)
allocate(tauatp(natmtot))
tauatp(:)=tau0atp
! initialise the previous total force on each atom
if (allocated(forcetotp)) deallocate(forcetotp)
allocate(forcetotp(3,natmtot))
forcetotp(:,:)=0.d0
! initial lattice vector step size
taulatv(:)=tau0latv
! initialise previous stress matrix
stressp(:,:)=0.d0
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
! open STRESSMAX.OUT and STRESS_OPT.OUT if required
if (latvopt.ne.0) then
  open(76,file='STRESSMAX.OUT',action='WRITE',form='FORMATTED')
  open(77,file='STRESS_OPT.OUT',action='WRITE',form='FORMATTED')
end if
if (mp_mpi) write(*,*)
do istp=1,maxlatvstp
  do jstp=1,maxatpstp
    if (mp_mpi) then
      write(*,'("Info(geomopt): atomic position optimisation step : ",I6)') jstp
    end if
! ground-state and forces calculation
    call gndstate
! subsequent calculations will read in the potential from STATE.OUT
    trdstate=.true.
! update the atomic positions
    call atpstep
! write total energy, forces, atomic positions, interatomic distances to file
    if (mp_mpi) then
      write(71,'(G22.12)') engytot
      call flushifc(71)
      write(72,'(G18.10)') forcemax
      call flushifc(72)
      write(73,*); write(73,*)
      write(73,'("! Lattice and atomic position optimisation steps : ",2I6)') &
       istp,jstp
      call writegeom(73)
      call flushifc(73)
      write(74,*); write(74,*)
      write(74,'("Lattice and atomic position optimisation steps : ",2I6)') &
       istp,jstp
      call writeiad(74)
      call flushifc(74)
      write(75,*); write(75,*)
      write(75,'("Lattice and atomic position optimisation steps : ",2I6)') &
       istp,jstp
      call writeforces(75)
      write(75,*)
      write(75,'("Maximum force magnitude over all atoms (target) : ",G18.10,&
       &" (",G18.10,")")') forcemax,epsforce
      call flushifc(75)
    end if
! check force convergence
    if (forcemax.le.epsforce) then
      if (mp_mpi) then
        write(75,*)
        write(75,'("Force convergence target achieved")')
      end if
      exit
    end if
    if (mp_mpi.and.(jstp.eq.maxatpstp)) then
      write(*,*)
      write(*,'("Warning(geomopt): atomic position optimisation failed to &
       &converge in ",I6," steps")') maxatpstp
    end if
! store the current forces array
    forcetotp(:,:)=forcetot(:,:)
! end loop over atomic position optimisation
  end do
! exit lattice optimisation loop if required
  if (latvopt.eq.0) exit
  if (mp_mpi) then
    write(*,'("Info(geomopt): lattice vector optimisation step : ",I6)') istp
  end if
! generate the stress matrix
  call genstress
! update the lattice vectors
  call latvstep
! write stress magnitude and matrix to file
  if (mp_mpi) then
    write(76,'(G18.10)') stressmax
    call flushifc(76)
    write(77,*)
    write(77,'("Lattice vector optimisation step : ",I6)') istp
    do j=1,3
      write(77,'(3G18.10)') (stress(i,j),i=1,3)
    end do
    call flushifc(77)
  end if
! check for stress convergence; stress may be non-zero because of volume
! constraint; checking change in stress matrix instead
  t1=sum(abs(stress(:,:)-stressp(:,:)))
  if (t1.le.epsstress*tau0latv) then
    if (mp_mpi) then
      write(77,*)
      write(77,'("Stress convergence target achieved")')
    end if
    exit
  end if
  if (mp_mpi.and.(istp.eq.maxlatvstp)) then
    write(*,*)
    write(*,'("Warning(geomopt): lattice vector optimisation failed to &
     &converge in ",I6," steps")') maxlatvstp
  end if
! store the current stress matrix
  stressp(:,:)=stress(:,:)
! end loop over lattice optimisation
end do
close(71); close(72); close(73); close(74); close(75)
if (latvopt.ne.0) then
  close(76); close(77)
end if
! ground-state should be run again after lattice vector optimisation
if (latvopt.ne.0) call gndstate
return
end subroutine

