
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvnlijji
! !INTERFACE:
subroutine genvnlijji(ikp,vnlijji)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced k-point set (in,integer)
!   vnlijji : non-local Coulomb matrix elements (out,real(nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Calculates the non-local Coulomb matrix elements of the type $(i-jj-i)$.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(out) :: vnlijji(nstsv,nstsv,nkpt)
! local variables
integer ngknr,ik,igk
integer ist1,ist2
integer ig,iq,igq0,iv(3)
real(8) cfq,v(3),t1
complex(8) zrho0,zt1
! allocatable arrays
integer, allocatable :: igkignr(:)
real(8), allocatable :: vgklnr(:,:),vgkcnr(:,:),gkcnr(:),tpgkcnr(:,:)
real(8), allocatable :: vgqc(:,:),tpgqc(:,:),gqc(:),jlgqr(:,:,:)
complex(8), allocatable :: sfacgknr(:,:),apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:)
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(igkignr(ngkmax))
allocate(vgklnr(3,ngkmax),vgkcnr(3,ngkmax),gkcnr(ngkmax),tpgkcnr(2,ngkmax))
allocate(vgqc(3,ngvec),tpgqc(2,ngvec),gqc(ngvec))
allocate(jlgqr(0:lnpsd+1,ngvec,nspecies))
allocate(sfacgknr(ngkmax,natmtot),apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxvr,ngvec),sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv),wfir2(ngrtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot),zvclir(ngrtot))
! factor for long-range term
cfq=0.5d0*(omega/pi)**2
! generate G+k-vectors
call gengpvec(vkl(:,ikp),vkc(:,ikp),ngknr,igkignr,vgklnr,vgkcnr)
! generate the spherical coordinates of the G+k-vectors
do igk=1,ngknr
  call sphcrd(vgkcnr(:,igk),gkcnr(igk),tpgkcnr(:,igk))
end do
! generate the structure factors
call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
! get the eigenvectors from file for non-reduced k-points
call getevecfv(vkl(:,ikp),vgklnr,evecfv)
call getevecsv(vkl(:,ikp),evecsv)
! calculate the wavefunctions for all states for passed non-reduced k-point ikp
call genwfsv(.false.,.false.,.false.,ngknr,igkignr,evalsv,apwalm,evecfv, &
 evecsv,wfmt2,ngrtot,wfir2)
! start loop over reduced k-point set
do ik=1,nkpt
! get the eigenvectors from file
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! calculate the wavefunctions for all states for the reduced k-point
  call genwfsv(.false.,.false.,.false.,ngk(1,ik),igkig(:,1,ik),evalsv,apwalm, &
   evecfv,evecsv,wfmt1,ngrtot,wfir1)
! determine q-vector
  iv(:)=ivk(:,ik)-ivk(:,ikp)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ik)-vkc(:,ikp)
  do ig=1,ngvec
! determine G+q-vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tpgqc(:,ig))
! spherical harmonics for G+q-vectors
    call genylm(lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
! compute the required spherical Bessel functions
  call genjlgpr(lnpsd+1,gqc,jlgqr)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
  do ist1=1,nstsv
    do ist2=1,nstsv
! calculate the complex overlap density
      call genzrho(.true.,wfmt2(:,:,:,:,ist2),wfmt1(:,:,:,:,ist1), &
       wfir2(:,:,ist2),wfir1(:,:,ist1),zrhomt,zrhoir)
! compute the potential and G=0 coefficient of the density
      call genzvclmt(nrcmt,nrcmtmax,rcmt,nrcmtmax,zrhomt,zvclmt)
      call zpotcoul(nrcmt,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq,sfacgq,zrhoir, &
       nrcmtmax,zvclmt,zvclir,zrho0)
      zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
      t1=cfq*wiq2(iq)*(dble(zrho0)**2+aimag(zrho0)**2)
      vnlijji(ist1,ist2,ik)=wkptnr*dble(zt1)+t1
! end loop over ist2
    end do
! end loop over ist1
  end do
! end loop over reduced k-point set
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr)
deallocate(vgqc,tpgqc,gqc,jlgqr)
deallocate(sfacgknr,apwalm,evecfv,evecsv,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine
!EOC

