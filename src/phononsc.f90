
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phononsc
use modmain
use modphonon
implicit none
! local variables
integer is,ia,ja,ias,jas
integer ip,nph,i,m,task0
real(8) dph,a,b,t1
real(8) ftp(3,maxatoms*maxspecies)
complex(8) zt1,zt2
complex(8) dyn(3,maxatoms,maxspecies)
! allocatable arrays
real(8), allocatable :: veffmtp(:,:,:)
real(8), allocatable :: veffirp(:)
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
if (allocated(dveffmt)) deallocate(dveffmt)
allocate(dveffmt(lmmaxvr,nrcmtmax,natmtot))
if (allocated(dveffir)) deallocate(dveffir)
allocate(dveffir(ngrtot))
! read original effective potential from file and make a copy
call readstate
if (allocated(veffmt0)) deallocate(veffmt0)
allocate(veffmt0(lmmaxvr,nrmtmax,natmtot))
if (allocated(veffir0)) deallocate(veffir0)
allocate(veffir0(ngrtot))
veffmt0(:,:,:)=veffmt(:,:,:)
veffir0(:)=veffir(:)
! store original parameters
natoms0(1:nspecies)=natoms(1:nspecies)
natmtot0=natmtot
avec0(:,:)=avec(:,:)
ainv0(:,:)=ainv(:,:)
atposc0(:,:,:)=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    atposc0(:,ia,is)=atposc(:,ia,is)
  end do
end do
ngrid0(:)=ngrid(:)
ngrtot0=ngrtot
task0=task
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
do m=0,nph
! restore input values
  natoms(1:nspecies)=natoms0(1:nspecies)
  avec(:,:)=avec0(:,:)
  atposc(:,:,:)=atposc0(:,:,:)
! generate the supercell
  call genphsc(m,deltaph)
! run the ground-state calculation
  call gndstate
! subsequent calculations will read in this potential
  task=1
! store the total force for the first displacement
  do ias=1,natmtot
    ftp(:,ias)=forcetot(:,ias)
  end do
! store the effective potential for the first displacement
  allocate(veffmtp(lmmaxvr,nrmtmax,natmtot))
  allocate(veffirp(ngrtot))
  veffmtp(:,:,:)=veffmt(:,:,:)
  veffirp(:)=veffir(:)
! restore input values
  natoms(1:nspecies)=natoms0(1:nspecies)
  avec(:,:)=avec0(:,:)
  atposc(:,:,:)=atposc0(:,:,:)
! generate the supercell again with twice the displacement
  dph=deltaph+deltaph
  call genphsc(m,dph)
! run the ground-state calculation again
  call gndstate
! compute the complex perturbing effective potential with implicit q-phase
  call phscdveff(m,veffmtp,veffirp)
  deallocate(veffmtp,veffirp)
! Fourier transform the force differences to obtain the dynamical matrix
  zt1=1.d0/(dble(nphsc)*deltaph)
! multiply by i for sin-like displacement
  if (m.eq.1) zt1=zt1*zi
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
          t1=-(forcetot(ip,jas)-ftp(ip,jas))
          dyn(ip,ia,is)=dyn(ip,ia,is)+zt2*t1
        end do
      end do
    end do
  end do
end do
! restore task number
task=task0
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
call writedveff
! delete the non-essential files
call phdelete
goto 10
end subroutine

