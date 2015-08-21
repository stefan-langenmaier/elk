
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengyk(ikp)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
! local variables
integer jkp,ik,jk,ist,jst
integer is,ia,ias,nrc,nrci
integer iv(3),ig,iq,igq0,m
real(8) cfq,v(3),t1
complex(8) zrho0,z1
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:),tpgqc(:,:),jlgqr(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
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
allocate(vgqc(3,ngvec),gqc(ngvec),tpgqc(2,ngvec))
allocate(jlgqr(0:lnpsd,ngvec,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxvr,ngvec),sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngtot,nspinor,nstsv),wfir2(ngtot,nspinor,nstsv))
allocate(wfcr(lmmaxvr,nrcmtmax,2),zfmt(lmmaxvr,nrcmtmax))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot),zvclir(ngtot))
! coefficient for long-range term
cfq=0.5d0*(omega/pi)**2
! get the eigenvectors from file for input k-point
call getevecfv(vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! equivalent reduced input k-point
jkp=ikmap(ivk(1,ikp),ivk(2,ikp),ivk(3,ikp))
! calculate the wavefunctions for occupied states of the input k-point
call genwfsv(.false.,.false.,.true.,ngk(1,ikp),igkig(:,1,ikp),occsv(:,jkp), &
 apwalm,evecfv,evecsv,wfmt1,ngtot,wfir1)
! start loop over non-reduced k-point set
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ikmap(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! determine q-vector
  iv(:)=ivk(:,ikp)-ivk(:,ik)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ikp)-vkc(:,ik)
  do ig=1,ngvec
! determine the G+q-vectors
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
  call genjlgpr(lnpsd,gqc,jlgqr)
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! get the eigenvectors from file for non-reduced k-points
  call getevecfv(vkl(:,ik),vgkl(:,:,1,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! calculate the wavefunctions for occupied states
  call genwfsv(.false.,.false.,.true.,ngk(1,ik),igkig(:,1,ik),occsv(:,jk), &
   apwalm,evecfv,evecsv,wfmt2,ngtot,wfir2)
!--------------------------------------------!
!    valence-valence-valence contribution    !
!--------------------------------------------!
  do jst=1,nstsv
    if (evalsv(jst,jk).lt.efermi) then
      do ist=1,nstsv
        if (evalsv(ist,jkp).lt.efermi) then
! calculate the complex overlap density
          call genzrho(.true.,spinpol,wfmt2(:,:,:,:,jst),wfmt1(:,:,:,:,ist), &
           wfir2(:,:,jst),wfir1(:,:,ist),zrhomt,zrhoir)
! calculate the Coulomb potential
          call genzvclmt(nrcmt,nrcmtmax,rcmt,nrcmtmax,zrhomt,zvclmt)
          call zpotcoul(nrcmt,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq,sfacgq, &
           zrhoir,nrcmtmax,zvclmt,zvclir,zrho0)
          z1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
          t1=cfq*wiq2(iq)*(dble(zrho0)**2+aimag(zrho0)**2)
!$OMP CRITICAL
          engyx=engyx-0.5d0*occmax*wkpt(ikp)*(wkptnr*dble(z1)+t1)
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
  nrci=nrcmtinr(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do jst=1,spnst(is)
      if (spcore(jst,is)) then
        do m=-spk(jst,is),spk(jst,is)-1
! pass m-1/2 to wavefcr
          call wavefcr(.false.,lradstp,is,ia,jst,m,nrcmtmax,wfcr)
          do ist=1,nstsv
            if (evalsv(ist,jkp).lt.efermi) then
! calculate the complex overlap density in spherical harmonics
              zfmt(:,1:nrc)=conjg(wfcr(:,1:nrc,1))*wfmt1(:,1:nrc,ias,1,ist)
              if (spinpol) then
                zfmt(:,1:nrc)=zfmt(:,1:nrc) &
                 +conjg(wfcr(:,1:nrc,2))*wfmt1(:,1:nrc,ias,2,ist)
              end if
              call zfsht(nrc,nrci,zfmt,zrhomt(:,:,ias))
! calculate the Coulomb potential
              call zpotclmt(lmaxvr,nrc,rcmt(:,is),lmmaxvr,zrhomt(:,:,ias), &
               zvclmt(:,:,ias))
              z1=zfmtinp(.true.,nrc,nrci,rcmt(:,is),zrhomt(:,:,ias), &
               zvclmt(:,:,ias))
!$OMP CRITICAL
              engyx=engyx-occmax*wkpt(ikp)*dble(z1)
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
deallocate(vgqc,gqc,tpgqc,jlgqr)
deallocate(evecfv,evecsv)
deallocate(apwalm,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2,wfcr)
deallocate(zrhomt,zrhoir,zvclmt,zvclir,zfmt)
return
end subroutine

