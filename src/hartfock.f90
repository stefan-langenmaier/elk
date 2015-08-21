
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hartfock
use modmain
use modmpi
implicit none
! local variables
logical exist
integer ik,is,ia,lp
real(8) etp,de
! allocatable arrays
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
call init2
! only the MPI master process should write files
if (mp_mpi) then
! open INFO.OUT file
  open(60,file='INFO'//trim(filext),action='WRITE',form='FORMATTED')
! open TOTENERGY.OUT
  open(61,file='TOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open FERMIDOS.OUT
  open(62,file='FERMIDOS'//trim(filext),action='WRITE',form='FORMATTED')
! open MOMENT.OUT if required
  if (spinpol) open(63,file='MOMENT'//trim(filext),action='WRITE', &
   form='FORMATTED')
! open FORCEMAX.OUT if required
  if (tforce) open(64,file='FORCEMAX'//trim(filext),action='WRITE', &
   form='FORMATTED')
! open DTOTENERGY.OUT
  open(66,file='DTOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! write out general information to INFO.OUT
  call writeinfo(60)
end if
! read the charge density and potentials from file
call readstate
! compute the effective potential
call poteff
! Fourier transform effective potential to G-space
call genveffig
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! generate the kinetic matrix elements in the Cartesian basis
call genkinmatc
! find the occupation numbers and Fermi energy
call occupy
10 continue
! set last self-consistent loop flag
tlast=.false.
etp=0.d0
! begin the self-consistent loop
if (mp_mpi) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
end if
do iscl=1,maxscl
  if (mp_mpi) then
    write(60,*)
    write(60,'("+--------------------+")')
    write(60,'("| Loop number : ",I4," |")') iscl
    write(60,'("+--------------------+")')
    call flushifc(60)
    write(*,*)
    write(*,'("Info(hartfock): self-consistent loop number : ",I4)') iscl
  end if
  if (iscl.ge.maxscl) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("Reached self-consistent loops maximum")')
    end if
    tlast=.true.
  end if
! compute the exchange-correlation potential for hybrid functionals
  if (hybmix.lt.1.d0) call potxc
! synchronise MPI processes
  call mpi_barrier(mpi_comm_world,ierror)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv)
!$OMP DO
  do ik=1,nkpt
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    allocate(evecsv(nstsv,nstsv))
    call getevecsv(vkl(:,ik),evecsv)
! solve the Hartree-Fock secular equation
    call seceqnhf(ik,evecsv)
! write the eigenvalues/vectors to file
    call putevalsv(ik,evalsv(:,ik))
    call putevecsv(ik,evecsv)
    deallocate(evecsv)
  end do
!$OMP END DO
!$OMP END PARALLEL
! synchronise MPI processes
  call mpi_barrier(mpi_comm_world,ierror)
! broadcast eigenvalue array to every process
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(evalsv(:,ik),nstsv,mpi_double_precision,lp,mpi_comm_world, &
     ierror)
  end do
! find the occupation numbers and Fermi energy
  call occupy
  if (mp_mpi) then
! write the occupation numbers to file
    do ik=1,nkpt
      call putoccsv(ik,occsv(:,ik))
    end do
! write out the eigenvalues and occupation numbers
    call writeeval
! write the Fermi energy to file
    call writefermi
  end if
! generate the density and magnetisation
  call rhomag
! compute the Coulomb potential
  call potcoul
! compute the energy components
  call energy
  if (mp_mpi) then
! output energy components
    call writeengy(60)
    write(60,*)
    write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
    write(60,'(" (states/Hartree/unit cell)")')
! write total energy to TOTENERGY.OUT and flush
    write(61,'(G22.12)') engytot
    call flushifc(61)
! write DOS at Fermi energy to FERMIDOS.OUT and flush
    write(62,'(G18.10)') fermidos
    call flushifc(62)
! output charges and moments
    call writechg(60)
! write total moment to MOMENT.OUT and flush
    if (spinpol) then
      write(63,'(3G18.10)') momtot(1:ndmag)
      call flushifc(63)
    end if
  end if
  if (tlast) goto 20
! compute the change in total energy and check for convergence
  if (iscl.ge.2) then
    de=abs(engytot-etp)
    if (mp_mpi) then
      write(60,*)
      write(60,'("Absolute change in total energy (target) : ",G18.10," (",&
       &G18.10,")")') de,epsengy
    end if
    if (de.lt.epsengy) then
      if (mp_mpi) then
        write(60,*)
        write(60,'("Energy convergence target achieved")')
      end if
      tlast=.true.
    end if
    if (mp_mpi) then
      write(66,'(G18.10)') de
      call flushifc(66)
    end if
  end if
  etp=engytot
! check for STOP file (only master process)
  if (mp_mpi) then
    inquire(file='STOP',exist=exist)
    if (exist) then
      write(60,*)
      write(60,'("STOP file exists - stopping self-consistent loop")')
      tlast=.true.
      open(50,file='STOP')
      close(50,status='DELETE')
    end if
  end if
! broadcast tlast from master process to all other processes
  call mpi_bcast(tlast,1,mpi_logical,0,mpi_comm_world,ierror)
end do
20 continue
if (mp_mpi) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop stopped |")')
  write(60,'("+------------------------------+")')
  if (maxscl.gt.1) then
    call writestate
    write(60,*)
    write(60,'("Wrote STATE.OUT")')
  end if
end if
!-----------------------!
!     compute forces    !
!-----------------------!
if (tforce) then
  call force
  if (mp_mpi) then
! output forces to INFO.OUT
    call writeforce(60)
! write maximum force magnitude to FORCEMAX.OUT
    write(64,'(G18.10)') forcemax
    call flushifc(64)
  end if
end if
!---------------------------------------!
!     perform structural relaxation     !
!---------------------------------------!
if (task.eq.6) then
  if (mp_mpi) then
    write(60,*)
    write(60,'("Maximum force magnitude (target) : ",G18.10," (",G18.10,")")') &
     forcemax,epsforce
    call flushifc(60)
  end if
! check force convergence
  if (forcemax.le.epsforce) then
    if (mp_mpi) then
      write(60,*)
      write(60,'("Force convergence target achieved")')
    end if
    goto 30
  end if
! update the atomic positions if forces are not converged
  call updatpos
  if (mp_mpi) then
! write optimised atomic positions and interatomic distances to file
    call writegeom(.true.)
    call writeiad(.true.)
    write(60,*)
    write(60,'("+--------------------------+")')
    write(60,'("| Updated atomic positions |")')
    write(60,'("+--------------------------+")')
    do is=1,nspecies
      write(60,*)
      write(60,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
      write(60,'(" atomic positions (lattice) :")')
      do ia=1,natoms(is)
        write(60,'(I4," : ",3F14.8)') ia,atposl(:,ia,is)
      end do
    end do
! add blank line to TOTENERGY.OUT, FERMIDOS.OUT, MOMENT.OUT and DTOTENERGY.OUT
    write(61,*)
    write(62,*)
    if (spinpol) write (63,*)
    write(66,*)
  end if
! begin new self-consistent loop with updated positions
  goto 10
end if
30 continue
if (mp_mpi) then
  write(60,*)
  write(60,'("+----------------------------+")')
  write(60,'("| Elk version ",I1.1,".",I1.1,".",I2.2," stopped |")') version
  write(60,'("+----------------------------+")')
! close the INFO.OUT file
  close(60)
! close the TOTENERGY.OUT file
  close(61)
! close the FERMIDOS.OUT file
  close(62)
! close the MOMENT.OUT file
  if (spinpol) close(63)
! close the FORCEMAX.OUT file
  if (tforce) close(64)
! close the DTOTENERGY.OUT file
  close(66)
end if
return
end subroutine

