
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengyk(ikp)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
! local variables
integer ik,ngknr,igk,ist,jst
integer is,ia,ias,nrc,m
integer iv(3),iq,ig,igq0
real(8) cfq,v(3),t1
complex(8) zrho0,zt1
! allocatable arrays
integer, allocatable :: igkignr(:)
real(8), allocatable :: vgklnr(:,:),vgkcnr(:,:),gkcnr(:),tpgkcnr(:,:)
real(8), allocatable :: vgqc(:,:),gqc(:),tpgqc(:,:),jlgqr(:,:,:)
real(8), allocatable :: evalsvp(:),evalsvnr(:)
complex(8), allocatable :: sfacgknr(:,:),apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: wfcr(:,:,:),zfmt(:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:)
! external functions
complex(8) zfinp,zfmtinp
external zfinp,zfmtinp
! allocate local arrays
allocate(igkignr(ngkmax))
allocate(vgklnr(3,ngkmax),vgkcnr(3,ngkmax),gkcnr(ngkmax),tpgkcnr(2,ngkmax))
allocate(vgqc(3,ngvec),gqc(ngvec),tpgqc(2,ngvec))
allocate(jlgqr(0:lnpsd+1,ngvec,nspecies))
allocate(evalsvp(nstsv),evalsvnr(nstsv))
allocate(sfacgknr(ngkmax,natmtot),apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxvr,ngvec),sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv),wfir2(ngrtot,nspinor,nstsv))
allocate(wfcr(lmmaxvr,nrcmtmax,2),zfmt(lmmaxvr,nrcmtmax))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot),zvclir(ngrtot))
! coefficient for long-range term
cfq=0.5d0*(omega/pi)**2
! get the eigenvalues/vectors from file for input k-point
call getevalsv(vkl(:,ikp),evalsvp)
call getevecfv(vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! calculate the wavefunctions for occupied states for the input k-point
call genwfsv(.false.,.false.,.true.,ngk(1,ikp),igkig(:,1,ikp),evalsvp,apwalm, &
 evecfv,evecsv,wfmt1,ngrtot,wfir1)
! start loop over non-reduced k-point set
do ik=1,nkptnr
! generate G+k-vectors
  call gengpvec(vkl(:,ik),vkc(:,ik),ngknr,igkignr,vgklnr,vgkcnr)
! generate the spherical coordinates of the G+k-vectors
  do igk=1,ngknr
    call sphcrd(vgkcnr(:,igk),gkcnr(igk),tpgkcnr(:,igk))
  end do
! generate the structure factors
  call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
  call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
! get the eigenvalues/vectors from file for non-reduced k-points
  call getevalsv(vkl(:,ik),evalsvnr)
  call getevecfv(vkl(:,ik),vgklnr,evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! determine q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ikp)-vkc(:,ik)
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
! calculate the wavefunctions for occupied states
  call genwfsv(.false.,.false.,.true.,ngknr,igkignr,evalsvnr,apwalm,evecfv, &
   evecsv,wfmt2,ngrtot,wfir2)
!--------------------------------------------!
!    valence-valence-valence contribution    !
!--------------------------------------------!
  do jst=1,nstsv
    if (evalsvnr(jst).lt.efermi) then
      do ist=1,nstsv
        if (evalsvp(ist).lt.efermi) then
! calculate the complex overlap density
          call genzrho(.true.,wfmt2(:,:,:,:,jst),wfmt1(:,:,:,:,ist), &
           wfir2(:,:,jst),wfir1(:,:,ist),zrhomt,zrhoir)
! calculate the Coulomb potential
          call genzvclmt(nrcmt,nrcmtmax,rcmt,nrcmtmax,zrhomt,zvclmt)
          call zpotcoul(nrcmt,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq,sfacgq, &
           zrhoir,nrcmtmax,zvclmt,zvclir,zrho0)
          zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
          t1=cfq*wiq2(iq)*(dble(zrho0)**2+aimag(zrho0)**2)
!$OMP CRITICAL
          engyx=engyx-0.5d0*occmax*wkpt(ikp)*(wkptnr*dble(zt1)+t1)
!$OMP END CRITICAL
! end loop over ist
        end if
      end do
! end loop over jst
    end if
  end do
! end loop over non-reduced k-point set
end do
!-----------------------------------------!
!    valence-core-valence contribution    !
!-----------------------------------------!
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do jst=1,spnst(is)
      if (spcore(jst,is)) then
        do m=-spk(jst,is),spk(jst,is)-1
! pass m-1/2 to wavefcr
          call wavefcr(lradstp,is,ia,jst,m,nrcmtmax,wfcr)
          do ist=1,nstsv
            if (evalsvp(ist).lt.efermi) then
! calculate the complex overlap density in spherical harmonics
              zfmt(:,1:nrc)=conjg(wfcr(:,1:nrc,1))*wfmt1(:,1:nrc,ias,1,ist)
              if (spinpol) then
                zfmt(:,1:nrc)=zfmt(:,1:nrc) &
                 +conjg(wfcr(:,1:nrc,2))*wfmt1(:,1:nrc,ias,2,ist)
              end if
              call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr, &
               zfmt,lmmaxvr,zzero,zrhomt(:,:,ias),lmmaxvr)
! calculate the Coulomb potential
              call zpotclmt(lmaxvr,nrc,rcmt(:,is),lmmaxvr,zrhomt(:,:,ias), &
               zvclmt(:,:,ias))
              zt1=zfmtinp(.true.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
               zrhomt(:,:,ias),zvclmt(:,:,ias))
!$OMP CRITICAL
              engyx=engyx-occmax*wkpt(ikp)*dble(zt1)
!$OMP END CRITICAL
! end loop over ist
            end if
          end do
! end loop over m
        end do
! end loop over jst
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr)
deallocate(vgqc,gqc,tpgqc,jlgqr)
deallocate(evalsvp,evalsvnr,evecfv,evecsv)
deallocate(sfacgknr,apwalm,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2,wfcr)
deallocate(zrhomt,zrhoir,zvclmt,zvclir,zfmt)
return
end subroutine

