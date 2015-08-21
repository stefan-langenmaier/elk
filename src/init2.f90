
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init2
use modmain
use modrdm
use modphonon
implicit none
! local variables
logical lsym(48)
integer isym,is,ia
integer ist,ic,m
real(8) ts0,ts1
real(8) boxl(3,4),t1

call timesec(ts0)

!---------------------!
!     q-point set     !
!---------------------!
! check if the system is an isolated molecule
if (molecule) ngridq(:)=1
! store the point group symmetries for reducing the q-point set
if (reduceq.eq.0) then
  nsymqpt=1
  symqpt(:,:,1)=symlat(:,:,1)
else
  lsym(:)=.false.
  do isym=1,nsymcrys
    lsym(lsplsymc(isym))=.true.
  end do
  nsymqpt=0
  do isym=1,nsymlat
    if (lsym(isym)) then
      nsymqpt=nsymqpt+1
      symqpt(:,:,nsymqpt)=symlat(:,:,isym)
    end if
  end do
end if
! OEP, Hartree-Fock, RPA epsilon, BSE or RDMFT
if ((xctype(1).lt.0).or.(task.eq.5).or.(task.eq.6).or.(task.eq.180).or. &
 (task.eq.185).or.(task.eq.300)) then
  ngridq(:)=ngridk(:)
end if
! allocate the q-point arrays
if (allocated(iqmap)) deallocate(iqmap)
allocate(iqmap(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
if (allocated(iqmapnr)) deallocate(iqmapnr)
allocate(iqmapnr(0:ngridq(1)-1,0:ngridq(2)-1,0:ngridq(3)-1))
nqptnr=ngridq(1)*ngridq(2)*ngridq(3)
if (allocated(ivq)) deallocate(ivq)
allocate(ivq(3,nqptnr))
if (allocated(vql)) deallocate(vql)
allocate(vql(3,nqptnr))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,nqptnr))
if (allocated(wqpt)) deallocate(wqpt)
allocate(wqpt(nqptnr))
! setup the q-point box (offset should always be zero)
boxl(:,:)=0.d0
boxl(1,2)=1.d0; boxl(2,3)=1.d0; boxl(3,4)=1.d0;
! generate the q-point set, note that the vectors vql and vqc are mapped to the
! first Brillouin zone
call genppts(.true.,nsymqpt,symqpt,ngridq,nqptnr,epslat,bvec,boxl,nqpt,iqmap, &
 iqmapnr,ivq,vql,vqc,wqpt,wqptnr)
! find the index for q = 0
do iq0=1,nqpt
  t1=sum(abs(vql(:,iq0)))
  if (t1.lt.epslat) goto 10
end do
write(*,*)
write(*,'("Error(init2): q = 0 not in q-point set")')
write(*,*)
stop
10 continue

!-----------------------------------------------!
!     OEP, Hartree-Fock and RDMFT variables     !
!-----------------------------------------------!
if ((xctype(1).lt.0).or.(task.eq.5).or.(task.eq.6).or.(task.eq.300)) then
! determine the 1/q^2 integral weights if required
  call genwiq2
! output the 1/q^2 integrals to WIQ2.OUT
  call writewiq2
end if
if (xctype(1).lt.0) then
! initialise OEP residual magnitude
  resoep=1.d0
! find maximum core states over all species
  ncrmax=0
  do is=1,nspecies
    do ia=1,natoms(is)
      ic=0
      do ist=1,spnst(is)
        if (spcore(ist,is)) then
          do m=-spk(ist,is),spk(ist,is)-1
            ic=ic+1
          end do
        end if
      end do
      ncrmax=max(ncrmax,ic)
    end do
  end do
! allocate and zero the complex exchange potential and field
  if (allocated(zvxmt)) deallocate(zvxmt)
  allocate(zvxmt(lmmaxvr,nrcmtmax,natmtot))
  zvxmt(:,:,:)=0.d0
  if (allocated(zvxir)) deallocate(zvxir)
  allocate(zvxir(ngrtot))
  zvxir(:)=0.d0
  if (spinpol) then
    if (allocated(zbxmt)) deallocate(zbxmt)
    allocate(zbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
    zbxmt(:,:,:,:)=0.d0
    if (allocated(zbxir)) deallocate(zbxir)
    allocate(zbxir(ngrtot,ndmag))
    zbxir(:,:)=0.d0
  end if
end if
if ((task.eq.5).or.(task.eq.6).or.(task.eq.300)) then
! allocate the kinetic matrix elements for Hartree-Fock/RDMFT
  if (allocated(kinmatc)) deallocate(kinmatc)
  allocate(kinmatc(nstsv,nstsv,nkpt))
end if
if (task.eq.300) then
  if (allocated(vclmat)) deallocate(vclmat)
  allocate(vclmat(nstsv,nstsv,nkpt))
  if (allocated(dkdc)) deallocate(dkdc)
  allocate(dkdc(nstsv,nstsv,nkpt))
end if
! Coulomb interaction regulator for RPA: 1/(q^2+lambda^2)
t1=twopi**3/omega
t1=t1/dble(nqptnr)
t1=t1**(1.d0/3.d0)
t1=t1/2.d0
lmda2rpa=t1**2

!--------------------------!
!     phonon variables     !
!--------------------------!
if (task.eq.202) then
  if (allocated(drhomt)) deallocate(drhomt)
  allocate(drhomt(lmmaxvr,nrcmtmax,natmtot))
  if (allocated(drhoir)) deallocate(drhoir)
  allocate(drhoir(ngrtot))
  if (allocated(dveffpw)) deallocate(dveffpw)
  allocate(dveffpw(ngrtot))
  if (allocated(dveffmt)) deallocate(dveffmt)
  allocate(dveffmt(lmmaxvr,nrcmtmax,natmtot))
  if (allocated(dveffir)) deallocate(dveffir)
  allocate(dveffir(ngrtot))
  if (allocated(dmagmt)) deallocate(dmagmt)
  if (allocated(dmagir)) deallocate(dmagir)
  if (allocated(dbxcpw)) deallocate(dbxcpw)
  if (allocated(dbxcmt)) deallocate(dbxcmt)
  if (allocated(dbxcir)) deallocate(dbxcir)
  if (spinpol) then
    allocate(dmagmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
    allocate(dmagir(ngrtot,ndmag))
    allocate(dbxcpw(ngrtot,ndmag))
    allocate(dbxcmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
    allocate(dbxcir(ngrtot,ndmag))
  end if
end if

call timesec(ts1)
timeinit=timeinit+ts1-ts0

return
end subroutine

