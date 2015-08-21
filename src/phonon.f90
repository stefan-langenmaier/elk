
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phonon
use modmain
use modphonon
implicit none
! local variables
logical tconv
integer is,ia,ja,ias,jas
integer ip,nph,i,p
real(8) dph,a,b,t1
real(8) forcetot1(3,maxatoms*maxspecies)
complex(8) zt1,zt2
complex(8) dyn(3,maxatoms,maxspecies)
! allocatable arrays
real(8), allocatable :: veffmt1(:,:,:),veffir1(:)
complex(8), allocatable :: dveffmt(:,:,:),dveffir(:)
!------------------------!
!     initialisation     !
!------------------------!
! require forces
tforce=.true.
! switch off automatic determination of muffin-tin radii
autormt=.false.
! no shifting of atomic basis allowed
tshift=.false.
! determine k-point grid size from radkpt
autokpt=.true.
! initialise universal variables
call init0
! initialise q-point dependent variables
call init2
! allocate the effective potential derivative arrays
allocate(dveffmt(lmmaxvr,nrcmtmax,natmtot),dveffir(ngrtot))
! store original parameters
natoms0(1:nspecies)=natoms(1:nspecies)
natmtot0=natmtot
avec0(:,:)=avec(:,:)
ainv0(:,:)=ainv(:,:)
binv0(:,:)=binv(:,:)
atposc0(:,:,:)=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    atposc0(:,ia,is)=atposc(:,ia,is)
  end do
end do
ngrid0(:)=ngrid(:)
ngrtot0=ngrtot
!---------------------------------------!
!     compute dynamical matrix rows     !
!---------------------------------------!
10 continue
natoms(1:nspecies)=natoms0(1:nspecies)
! find a dynamical matrix to calculate
call dyntask(80)
! if nothing more to do then reset input values and return
if (iqph.eq.0) then
  call readinput
  return
end if
write(*,'("Info(phonon): working on ",A)') 'DYN'//trim(filext)
! phonon dry run: just generate empty DYN files
if (task.eq.205) goto 10
dyn(:,:,:)=0.d0
dveffmt(:,:,:)=0.d0
dveffir(:)=0.d0
! check to see if mass is considered infinite
if (spmass(isph).le.0.d0) goto 20
! loop over phases: 0 = cos and 1 = sin displacements
if ((ivq(1,iqph).eq.0).and.(ivq(2,iqph).eq.0).and.(ivq(3,iqph).eq.0)) then
  nph=0
else
  nph=1
end if
! initial supercell density constructed from atomic densities
trdstate=.false.
! flag for checking if ground-state calculations converged
tconv=.true.
! loop over cos and sin displacements
do p=0,nph
! restore input values
  natoms(1:nspecies)=natoms0(1:nspecies)
  avec(:,:)=avec0(:,:)
  atposc(:,:,:)=atposc0(:,:,:)
! generate the supercell
  call genphsc(p,deltaph)
! run the ground-state calculation
  call gndstate
  if (iscl.ge.maxscl) tconv=.false.
! subsequent calculations will read in this supercell potential
  trdstate=.true.
! store the total force for the first displacement
  do ias=1,natmtot
    forcetot1(:,ias)=forcetot(:,ias)
  end do
! store the effective potential for the first displacement
  allocate(veffmt1(lmmaxvr,nrmtmax,natmtot),veffir1(ngrtot))
  veffmt1(:,:,:)=veffmt(:,:,:)
  veffir1(:)=veffir(:)
! restore input values
  natoms(1:nspecies)=natoms0(1:nspecies)
  avec(:,:)=avec0(:,:)
  atposc(:,:,:)=atposc0(:,:,:)
! generate the supercell again with twice the displacement
  dph=deltaph+deltaph
  call genphsc(p,dph)
! run the ground-state calculation again
  call gndstate
  if (iscl.ge.maxscl) tconv=.false.
! compute the complex effective potential derivative with implicit q-phase
  call phdveff(p,veffmt1,veffir1,dveffmt,dveffir)
  deallocate(veffmt1,veffir1)
! Fourier transform the force differences to obtain the dynamical matrix
  zt1=1.d0/(dble(nphsc)*deltaph)
! multiply by i for sin-like displacement
  if (p.eq.1) zt1=zt1*zi
  jas=0
  do is=1,nspecies
    ja=0
    do ia=1,natoms0(is)
      do i=1,nphsc
        ja=ja+1
        jas=jas+1
        t1=-dot_product(vqc(:,iqph),vphsc(:,i))
        zt2=zt1*cmplx(cos(t1),sin(t1),8)
        do ip=1,3
          t1=-(forcetot(ip,jas)-forcetot1(ip,jas))
          dyn(ip,ia,is)=dyn(ip,ia,is)+zt2*t1
        end do
      end do
    end do
  end do
end do
if (.not.tconv) then
  write(*,*)
  write(*,'("Warning(phonon): supercell calculation failed to converge")')
  write(*,'("Current dynamical matrix row probably inaccurate")')
end if
20 continue
! write dynamical matrix row to file
do is=1,nspecies
  do ia=1,natoms0(is)
    do ip=1,3
      a=dble(dyn(ip,ia,is))
      b=aimag(dyn(ip,ia,is))
      if (abs(a).lt.1.d-12) a=0.d0
      if (abs(b).lt.1.d-12) b=0.d0
      write(80,'(2G18.10," : is = ",I4,", ia = ",I4,", ip = ",I4)') a,b,is,ia,ip
    end do
  end do
end do
close(80)
! write the complex effective potential derivative to file
call writedveff(dveffmt,dveffir)
! delete the non-essential files
call phdelete
goto 10
end subroutine

