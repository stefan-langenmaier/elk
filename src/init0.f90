
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: init0
! !INTERFACE:
subroutine init0
! !USES:
use modmain
use modxcifc
use modldapu
use modtest
! !DESCRIPTION:
!   Performs basic consistency checks as well as allocating and initialising
!   global variables not dependent on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ist
integer l,m,lm
real(8) rsum,t1
real(8) ts0,ts1

!-------------------------------!
!     zero timing variables     !
!-------------------------------!
timeinit=0.d0
timemat=0.d0
timefv=0.d0
timesv=0.d0
timerho=0.d0
timepot=0.d0
timefor=0.d0
call timesec(ts0)

!------------------------------------!
!     angular momentum variables     !
!------------------------------------!
lmmaxvr=(lmaxvr+1)**2
lmmaxapw=(lmaxapw+1)**2
lmmaxmat=(lmaxmat+1)**2
lmmaxinr=(lmaxinr+1)**2
if (lmaxvr.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(init0): lmaxvr > lmaxapw : ",2I8)') lmaxvr,lmaxapw
  write(*,*)
  stop
end if
if (lmaxmat.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(init0): lmaxmat > lmaxapw : ",2I8)') lmaxmat,lmaxapw
  write(*,*)
  stop
end if
! check DOS lmax is within range
lmaxdos=min(lmaxdos,lmaxapw)
! index to (l,m) pairs
if (allocated(idxlm)) deallocate(idxlm)
allocate(idxlm(0:lmaxapw,-lmaxapw:lmaxapw))
if (allocated(idxil)) deallocate(idxil)
allocate(idxil(lmmaxapw))
if (allocated(idxim)) deallocate(idxim)
allocate(idxim(lmmaxapw))
lm=0
do l=0,lmaxapw
  do m=-l,l
    lm=lm+1
    idxlm(l,m)=lm
    idxil(lm)=l
    idxim(lm)=m
  end do
end do
! array of i^l and (-i)^l values
if (allocated(zil)) deallocate(zil)
if (allocated(zilc)) deallocate(zilc)
allocate(zil(0:lmaxapw),zilc(0:lmaxapw))
do l=0,lmaxapw
  zil(l)=zi**l
  zilc(l)=conjg(zil(l))
end do

!------------------------------------!
!     index to atoms and species     !
!------------------------------------!
natmmax=0
ias=0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=ias+1
    idxas(ia,is)=ias
    idxis(ias)=is
    idxia(ias)=ia
  end do
! maximum number of atoms over all species
  natmmax=max(natmmax,natoms(is))
end do
! total number of atoms
natmtot=ias

!------------------------!
!     spin variables     !
!------------------------!
if (spinsprl) then
  spinpol=.true.
  spinorb=.false.
  select case(task)
  case(51,52,53,61,62,63)
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with task ",I4)') task
    write(*,*)
    stop
  end select
  if (xctype(1).lt.0) then
    write(*,*)
    write(*,'("Error(init0): spin-spirals do not work with the OEP method")')
    write(*,*)
    stop
  end if
end if
! spin-orbit coupling or fixed spin moment implies spin-polarised calculation
if ((spinorb).or.(fixspin.ne.0).or.(spinsprl)) spinpol=.true.
! number of spinor components and maximum allowed occupancy
if (spinpol) then
  nspinor=2
  occmax=1.d0
else
  nspinor=1
  occmax=2.d0
end if
! number of spin-dependent first-variational functions per state and map from
! second- to first-variational spin index
if (spinsprl) then
  nspnfv=2
  jspnfv(1)=1
  jspnfv(2)=2
else
  nspnfv=1
  jspnfv(1)=1
  jspnfv(2)=1
end if
! spin-polarised calculations require second-variational eigenvectors
if (spinpol) tevecsv=.true.
! Hartree-Fock/RDMFT requires second-variational eigenvectors
if ((task.eq.5).or.(task.eq.300)) tevecsv=.true.
! get exchange-correlation functional data
call getxcdata(xctype,xcdescr,xcspin,xcgrad,hybrid,hybridc)
if ((spinpol).and.(xcspin.eq.0)) then
  write(*,*)
  write(*,'("Error(init0): requested spin-polarised run with &
   &spin-unpolarised")')
  write(*,'(" exchange-correlation functional")')
  write(*,*)
  stop
end if
! check for collinearity in the z-direction and set the dimension of the
! magnetisation and exchange-correlation vector fields
if (spinpol) then
  ndmag=1
  if ((abs(bfieldc0(1)).gt.epslat).or.(abs(bfieldc0(2)).gt.epslat)) ndmag=3
  do is=1,nspecies
    do ia=1,natoms(is)
      if ((abs(bfcmt0(1,ia,is)).gt.epslat).or. &
          (abs(bfcmt0(2,ia,is)).gt.epslat)) ndmag=3
    end do
  end do
! source-free fields and spin-spirals are non-collinear in general
  if ((nosource).or.(spinsprl)) ndmag=3
! spin-orbit coupling is non-collinear in general
  if (spinorb) ndmag=3
else
  ndmag=0
end if
! set the non-collinear flag
if (ndmag.eq.3) then
  ncmag=.true.
else
  ncmag=.false.
end if
! check for meta-GGA with non-collinearity
if ((xcgrad.eq.3).and.ncmag) then
  write(*,*)
  write(*,'("Error(init0): meta-GGA is not valid for non-collinear magnetism")')
  write(*,*)
  stop
end if
! spin-polarised cores
if (.not.spinpol) spincore=.false.
if (fixspin.ne.0) then
! set fixed spin moment effective field to zero
  bfsmc(:)=0.d0
! set muffin-tin FSM fields to zero
  if (allocated(bfsmcmt)) deallocate(bfsmcmt)
  allocate(bfsmcmt(3,natmtot))
  bfsmcmt(:,:)=0.d0
end if
! number of independent spin components of the f_xc spin tensor
if (spinpol) then
  if (ncmag) then
    nscfxc=10
  else
    nscfxc=3
  end if
else
  nscfxc=1
end if
! set the magnetic fields to the initial values
bfieldc(:)=bfieldc0(:)
bfcmt(:,:,:)=bfcmt0(:,:,:)
! if reducebf < 1 then reduce the external magnetic fields immediately for
! non-self-consistent calculations or resumptions
if (reducebf.lt.1.d0-epslat) then
  if ((task.ge.10).and.(task.ne.28).and.(task.ne.200).and.(task.ne.201).and. &
   (task.ne.350).and.(task.ne.351)) then
    bfieldc(:)=0.d0
    bfcmt(:,:,:)=0.d0
  end if
end if

!----------------------------------!
!     crystal structure set up     !
!----------------------------------!
! generate the reciprocal lattice vectors and unit cell volume
call reciplat
! inverse of the lattice vector matrix
call r3minv(avec,ainv)
! inverse of the reciprocal vector matrix
call r3minv(bvec,binv)
! Cartesian coordinates of the spin-spiral vector
call r3mv(bvec,vqlss,vqcss)
do is=1,nspecies
  do ia=1,natoms(is)
! map atomic lattice coordinates to [0,1)
    call r3frac(epslat,atposl(:,ia,is))
! determine atomic Cartesian coordinates
    call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
  end do
end do
! check muffin-tins are not too close together
call checkmt

!-------------------------------!
!     vector fields E and A     !
!-------------------------------!
efieldpol=.false.
if (sum(abs(efieldc(:))).gt.epslat) then
  efieldpol=.true.
  tshift=.false.
! electric field vector in lattice coordinates
  call r3mv(ainv,efieldc,efieldl)
end if
afieldpol=.false.
if (sum(abs(afieldc(:))).gt.epslat) then
  afieldpol=.true.
! vector potential added in second-variational step
  tevecsv=.true.
end if

!---------------------------------!
!     crystal symmetry set up     !
!---------------------------------!
call symmetry

!-----------------------!
!     radial meshes     !
!-----------------------!
nrmtmax=1
nrcmtmax=1
do is=1,nspecies
! make the muffin-tin mesh commensurate with lradstp
  nrmt(is)=nrmt(is)-mod(nrmt(is)-1,lradstp)
  nrmtmax=max(nrmtmax,nrmt(is))
! number of coarse radial mesh points
  nrcmt(is)=(nrmt(is)-1)/lradstp+1
  nrcmtmax=max(nrcmtmax,nrcmt(is))
end do
! set up atomic and muffin-tin radial meshes
call genrmesh

!--------------------------------------!
!     charges and number of states     !
!--------------------------------------!
chgzn=0.d0
chgcrtot=0.d0
chgval=0.d0
spnstmax=0
nstcr=0
do is=1,nspecies
! nuclear charge
  chgzn=chgzn+spzn(is)*natoms(is)
! find the maximum number of atomic states
  spnstmax=max(spnstmax,spnst(is))
! compute the electronic charge for each species, as well as the total core and
! valence charge
  spze(is)=0.d0
  chgcr(is)=0.d0
  do ist=1,spnst(is)
    spze(is)=spze(is)+spocc(ist,is)
    if (spcore(ist,is)) then
      chgcr(is)=chgcr(is)+spocc(ist,is)
      nstcr=nstcr+2*spk(ist,is)*natoms(is)
    else
      chgval=chgval+spocc(ist,is)*natoms(is)
    end if
  end do
  chgcrtot=chgcrtot+chgcr(is)*natoms(is)
end do
! add excess charge
chgval=chgval+chgexs
! total charge
chgtot=chgcrtot+chgval
if (chgtot.lt.1.d-8) then
  write(*,*)
  write(*,'("Error(init0): zero total charge")')
  write(*,*)
  stop
end if
! effective Wigner radius
rwigner=(3.d0/(fourpi*(chgtot/omega)))**(1.d0/3.d0)

!-------------------------!
!     G-vector arrays     !
!-------------------------!
if (nspecies.eq.0) isgkmax=-2
! determine gkmax from rgkmax and the muffin-tin radius
if (isgkmax.eq.-2) then
  gkmax=rgkmax/2.d0
else
  if ((isgkmax.ge.1).and.(isgkmax.le.nspecies)) then
! use user-specified muffin-tin radius
    gkmax=rgkmax/rmt(isgkmax)
  else if (isgkmax.eq.-1) then
! use average muffin-tin radius
    rsum=0.d0
    do is=1,nspecies
      rsum=rsum+dble(natoms(is))*rmt(is)
    end do
    rsum=rsum/dble(natmtot)
    gkmax=rgkmax/rsum
  else
! use minimum muffin-tin radius (isgkmax=-3)
    gkmax=rgkmax/minval(rmt(1:nspecies))
  end if
end if
! ensure |G| cut-off is at least twice |G+k| cut-off
gmaxvr=max(gmaxvr,2.d0*gkmax+epslat)
! find the G-vector grid sizes
call gridsize(avec,gmaxvr,ngridg,ngtot,intgv)
! allocate global G-vector arrays
if (allocated(ivg)) deallocate(ivg)
allocate(ivg(3,ngtot))
if (allocated(ivgig)) deallocate(ivgig)
allocate(ivgig(intgv(1,1):intgv(2,1),intgv(1,2):intgv(2,2), &
 intgv(1,3):intgv(2,3)))
if (allocated(igfft)) deallocate(igfft)
allocate(igfft(ngtot))
if (allocated(vgc)) deallocate(vgc)
allocate(vgc(3,ngtot))
if (allocated(gc)) deallocate(gc)
allocate(gc(ngtot))
! generate the G-vectors
call gengvec(ngridg,ngtot,intgv,bvec,gmaxvr,ngvec,ivg,ivgig,igfft,vgc,gc)
! write number of G-vectors to test file
call writetest(900,'number of G-vectors',iv=ngvec)
! Poisson solver pseudocharge density constant
if (nspecies.gt.0) then
  t1=0.25d0*gmaxvr*maxval(rmt(1:nspecies))
else
  t1=0.25d0*gmaxvr*2.d0
end if
npsd=max(nint(t1),1)
lnpsd=lmaxvr+npsd+1
! compute the spherical Bessel functions j_l(|G|R_mt)
if (allocated(jlgr)) deallocate(jlgr)
allocate(jlgr(0:lnpsd,ngvec,nspecies))
call genjlgpr(lnpsd,gc,jlgr)
! generate the spherical harmonics of the G-vectors
call genylmg
! allocate structure factor array for G-vectors
if (allocated(sfacg)) deallocate(sfacg)
allocate(sfacg(ngvec,natmtot))
! generate structure factors for G-vectors
call gensfacgp(ngvec,vgc,ngvec,sfacg)
! generate the smooth step function form factors
if (allocated(ffacg)) deallocate(ffacg)
allocate(ffacg(ngtot,nspecies))
do is=1,nspecies
  call genffacgp(is,gc,ffacg(:,is))
end do
! generate the characteristic function
call gencfun

!-------------------------!
!     atoms and cores     !
!-------------------------!
! solve the Kohn-Sham-Dirac equations for all atoms
call allatoms
! allocate core state occupancy and eigenvalue arrays and set to default
if (allocated(occcr)) deallocate(occcr)
allocate(occcr(spnstmax,natmtot))
if (allocated(evalcr)) deallocate(evalcr)
allocate(evalcr(spnstmax,natmtot))
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist=1,spnst(is)
      occcr(ist,ias)=spocc(ist,is)
      evalcr(ist,ias)=speval(ist,is)
    end do
  end do
end do
! allocate core state radial wavefunction array
if (allocated(rwfcr)) deallocate(rwfcr)
allocate(rwfcr(spnrmax,2,spnstmax,natmtot))
! number of core spin channels
if (spincore) then
  nspncr=2
else
  nspncr=1
end if
! allocate core state charge density array
if (allocated(rhocr)) deallocate(rhocr)
allocate(rhocr(spnrmax,natmtot,nspncr))

!---------------------------------------!
!     charge density and potentials     !
!---------------------------------------!
! allocate charge density arrays
if (allocated(rhomt)) deallocate(rhomt)
allocate(rhomt(lmmaxvr,nrmtmax,natmtot))
if (allocated(rhoir)) deallocate(rhoir)
allocate(rhoir(ngtot))
! allocate magnetisation arrays
if (allocated(magmt)) deallocate(magmt)
if (allocated(magir)) deallocate(magir)
if (spinpol) then
  allocate(magmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(magir(ngtot,ndmag))
end if
! Coulomb potential
if (allocated(vclmt)) deallocate(vclmt)
allocate(vclmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(vclir)) deallocate(vclir)
allocate(vclir(ngtot))
! exchange energy density
if (allocated(exmt)) deallocate(exmt)
allocate(exmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(exir)) deallocate(exir)
allocate(exir(ngtot))
! correlation energy density
if (allocated(ecmt)) deallocate(ecmt)
allocate(ecmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(ecir)) deallocate(ecir)
allocate(ecir(ngtot))
! exchange-correlation potential
if (allocated(vxcmt)) deallocate(vxcmt)
allocate(vxcmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(vxcir)) deallocate(vxcir)
allocate(vxcir(ngtot))
! effective Kohn-Sham potential
if (allocated(vsmt)) deallocate(vsmt)
allocate(vsmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(vsir)) deallocate(vsir)
allocate(vsir(ngtot))
if (allocated(vsig)) deallocate(vsig)
allocate(vsig(ngvec))
! exchange-correlation magnetic and Kohn-Sham effective fields
if (allocated(bxcmt)) deallocate(bxcmt)
if (allocated(bxcir)) deallocate(bxcir)
if (allocated(bsmt)) deallocate(bsmt)
if (allocated(bsir)) deallocate(bsir)
if (spinpol) then
  allocate(bxcmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(bxcir(ngtot,ndmag))
  allocate(bsmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
  allocate(bsir(ngtot,ndmag))
end if
! spin-orbit coupling radial function
if (allocated(socfr)) deallocate(socfr)
if (spinorb) then
  allocate(socfr(nrcmtmax,natmtot))
end if
! allocate muffin-tin charge and moment arrays
if (allocated(chgcrlk)) deallocate(chgcrlk)
allocate(chgcrlk(natmtot))
if (allocated(chgmt)) deallocate(chgmt)
allocate(chgmt(natmtot))
if (allocated(mommt)) deallocate(mommt)
allocate(mommt(3,natmtot))

!-------------------------!
!     force variables     !
!-------------------------!
if (allocated(forcehf)) deallocate(forcehf)
allocate(forcehf(3,natmtot))
if (allocated(forceibs)) deallocate(forceibs)
allocate(forceibs(3,natmtot))
if (allocated(forcetot)) deallocate(forcetot)
allocate(forcetot(3,natmtot))

!-------------------------!
!     LDA+U variables     !
!-------------------------!
if ((ldapu.ne.0).or.(task.eq.17)) then
! LDA+U requires second-variational eigenvectors
  tevecsv=.true.
! density matrices
  if (allocated(dmatlu)) deallocate(dmatlu)
  allocate(dmatlu(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
! potential matrix elements
  if (allocated(vmatlu)) deallocate(vmatlu)
  allocate(vmatlu(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
! zero the potential
  vmatlu(:,:,:,:,:)=0.d0
! energy for each atom
  if (allocated(engyalu)) deallocate(engyalu)
  allocate(engyalu(natmtot))
! interpolation constants (alpha)
  if (allocated(alphalu)) deallocate(alphalu)
  allocate(alphalu(natmtot))
end if

!-----------------------!
!     miscellaneous     !
!-----------------------!
! determine nuclear radii and volumes
call nuclei
! determine the nuclear-nuclear energy
call energynn
! get smearing function description
call getsdata(stype,sdescr)
! get mixing type description
call getmixdata(mixtype,mixdescr)
! generate the spherical harmonic transform (SHT) matrices
call genshtmat
! allocate 1D plotting arrays
if (allocated(dvp1d)) deallocate(dvp1d)
allocate(dvp1d(nvp1d))
if (allocated(vplp1d)) deallocate(vplp1d)
allocate(vplp1d(3,npp1d))
if (allocated(dpp1d)) deallocate(dpp1d)
allocate(dpp1d(npp1d))
! zero self-consistent loop number
iscl=0
tlast=.false.
! set the Fermi energy to zero
efermi=0.d0
! input q-vector in Cartesian coordinates
call r3mv(bvec,vecql,vecqc)

call timesec(ts1)
timeinit=timeinit+ts1-ts0

return
end subroutine
!EOC

