
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepvnlk(ikp,vnlcv,vnlvv)
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(out) :: vnlcv(ncrmax,natmtot,nstsv)
complex(8), intent(out) :: vnlvv(nstsv,nstsv)
! local variables
integer ngknr,ik,jk,igk
integer ist1,ist2,ist3
integer is,ia,ias,nrc
integer ic,jc,m1,m2
integer iq,ig,iv(3),igq0
real(8) v(3),cfq
complex(8) zrho01,zrho02,zt1,zt2
! automatic arrays
complex(8) sfacgq0(natmtot)
! allocatable arrays
integer, allocatable :: igkignr(:)
real(8), allocatable :: vgklnr(:,:),vgkcnr(:,:),gkcnr(:),tpgkcnr(:,:)
real(8), allocatable :: vgqc(:,:),tpgqc(:,:),gqc(:)
real(8), allocatable :: jlgqr(:,:,:),jlgq0r(:,:,:)
complex(8), allocatable :: sfacgknr(:,:),apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: wfcr1(:,:,:),wfcr2(:,:,:)
complex(8), allocatable :: zrhomt1(:,:,:,:),zrhomt2(:,:,:),zrhoir1(:,:)
complex(8), allocatable :: zvclmt(:,:,:),zvclir(:),zfmt(:,:)
! external functions
complex(8) zfinp,zfmtinp
external zfinp,zfmtinp
! allocate local arrays
allocate(igkignr(ngkmax))
allocate(vgklnr(3,ngkmax),vgkcnr(3,ngkmax),gkcnr(ngkmax),tpgkcnr(2,ngkmax))
allocate(vgqc(3,ngvec),tpgqc(2,ngvec),gqc(ngvec))
allocate(jlgqr(0:lnpsd+1,ngvec,nspecies),jlgq0r(0:lmaxvr,nrcmtmax,nspecies))
allocate(sfacgknr(ngkmax,natmtot),apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxvr,ngvec),sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv),wfir2(ngrtot,nspinor,nstsv))
allocate(wfcr1(lmmaxvr,nrcmtmax,2),wfcr2(lmmaxvr,nrcmtmax,2))
allocate(zrhomt1(lmmaxvr,nrcmtmax,natmtot,nstsv))
allocate(zrhomt2(lmmaxvr,nrcmtmax,nstcr))
allocate(zrhoir1(ngrtot,nstsv))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot),zvclir(ngrtot))
allocate(zfmt(lmmaxvr,nrcmtmax))
! factor for long-range term
cfq=0.5d0*(omega/pi)**2
! zero the Coulomb matrix elements
vnlcv(:,:,:)=0.d0
vnlvv(:,:)=0.d0
! get the eigenvectors from file for input k-point
call getevecfv(vkl(:,ikp),vgkl(:,:,:,ikp),evecfv)
call getevecsv(vkl(:,ikp),evecsv)
! find the matching coefficients
call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! calculate the wavefunctions for all states for the input k-point
call genwfsv(.false.,.false.,.false.,ngk(1,ikp),igkig(:,1,ikp),evalsv,apwalm, &
 evecfv,evecsv,wfmt1,ngrtot,wfir1)
! start loop over non-reduced k-point set
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ikmap(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! generate G+k-vectors
  call gengpvec(vkl(:,ik),vkc(:,ik),ngknr,igkignr,vgklnr,vgkcnr)
! generate the spherical coordinates of the G+k-vectors
  do igk=1,ngknr
    call sphcrd(vgkcnr(:,igk),gkcnr(igk),tpgkcnr(:,igk))
  end do
! get the eigenvectors from file for non-reduced k-points
  call getevecfv(vkl(:,ik),vgklnr,evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! generate the structure factors
  call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
  call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
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
! spherical harmonics for G+q-vector
    call genylm(lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
  sfacgq0(:)=sfacgq(igq0,:)
! compute the required spherical Bessel functions
  call genjlgpr(lnpsd+1,gqc,jlgqr)
  call genjlgq0r(gqc(igq0),jlgq0r)
! calculate the wavefunctions for occupied states
  call genwfsv(.false.,.false.,.true.,ngknr,igkignr,evalsv(:,jk),apwalm, &
   evecfv,evecsv,wfmt2,ngrtot,wfir2)
  do ist3=1,nstsv
    if (evalsv(ist3,jk).lt.efermi) then
! compute the complex overlap densities for all valence-valence states
      do ist1=1,nstsv
        call genzrho(.true.,wfmt2(:,:,:,:,ist3),wfmt1(:,:,:,:,ist1), &
         wfir2(:,:,ist3),wfir1(:,:,ist1),zrhomt1(:,:,:,ist1),zrhoir1(:,ist1))
      end do
! compute the complex overlap densities for all valence-core states
      jc=0
      do is=1,nspecies
        nrc=nrcmt(is)
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          do ist1=1,spnst(is)
            if (spcore(ist1,is)) then
              do m1=-spk(ist1,is),spk(ist1,is)-1
                jc=jc+1
! pass m-1/2 to wavefcr
                call wavefcr(lradstp,is,ia,ist1,m1,nrcmtmax,wfcr1)
                zfmt(:,1:nrc)=conjg(wfmt2(:,1:nrc,ias,1,ist3))*wfcr1(:,1:nrc,1)
                if (spinpol) then
                  zfmt(:,1:nrc)=zfmt(:,1:nrc)+conjg(wfmt2(:,1:nrc,ias,2,ist3)) &
                   *wfcr1(:,1:nrc,2)
                end if
                call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr, &
                 zfmt,lmmaxvr,zzero,zrhomt2(:,:,jc),lmmaxvr)
              end do
            end if
          end do
        end do
      end do
      do ist2=1,nstsv
        if (evalsv(ist2,ikp).gt.efermi) then
! calculate the Coulomb potential
          call genzvclmt(nrcmt,nrcmtmax,rcmt,nrcmtmax,zrhomt1(:,:,:,ist2), &
           zvclmt)
          call zpotcoul(nrcmt,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq,sfacgq, &
           zrhoir1(:,ist2),nrcmtmax,zvclmt,zvclir,zrho02)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
          do ist1=1,nstsv
            if (evalsv(ist1,ikp).lt.efermi) then
              zt1=zfinp(.true.,zrhomt1(:,:,:,ist1),zvclmt,zrhoir1(:,ist1), &
               zvclir)
! compute the density coefficient of the smallest G+q-vector
              call zrhogp(jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt1(:,:,:,ist1), &
               zrhoir1(:,ist1),zrho01)
              zt2=cfq*wiq2(iq)*(conjg(zrho01)*zrho02)
              vnlvv(ist1,ist2)=vnlvv(ist1,ist2)-(wkptnr*zt1+zt2)
            end if
          end do
!-------------------------------------------!
!     core-valence-valence contribution     !
!-------------------------------------------!
          jc=0
          do is=1,nspecies
            nrc=nrcmt(is)
            do ia=1,natoms(is)
              ias=idxas(ia,is)
              ic=0
              do ist1=1,spnst(is)
                if (spcore(ist1,is)) then
                  do m1=-spk(ist1,is),spk(ist1,is)-1
                    ic=ic+1
                    jc=jc+1
                    zt1=zfmtinp(.true.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
                     zrhomt2(:,:,jc),zvclmt(:,:,ias))
                    vnlcv(ic,ias,ist2)=vnlcv(ic,ias,ist2)-wkptnr*zt1
                  end do
! end loop over ist1
                end if
              end do
! end loops over atoms and species
            end do
          end do
! end loop over ist2
        end if
      end do
! end loop over ist3
    end if
  end do
! end loop over non-reduced k-point set
end do
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist3=1,spnst(is)
      if (spcore(ist3,is)) then
        do m1=-spk(ist3,is),spk(ist3,is)-1
! pass m-1/2 to wavefcr
          call wavefcr(lradstp,is,ia,ist3,m1,nrcmtmax,wfcr1)
! compute the complex overlap densities for the core-valence states
          do ist1=1,nstsv
            zfmt(:,1:nrc)=conjg(wfcr1(:,1:nrc,1))*wfmt1(:,1:nrc,ias,1,ist1)
            if (spinpol) then
              zfmt(:,1:nrc)=zfmt(:,1:nrc)+conjg(wfcr1(:,1:nrc,2)) &
               *wfmt1(:,1:nrc,ias,2,ist1)
            end if
            call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr, &
             zfmt,lmmaxvr,zzero,zrhomt1(:,:,ias,ist1),lmmaxvr)
          end do
! compute the complex overlap densities for the core-core states
          ic=0
          do ist1=1,spnst(is)
            if (spcore(ist1,is)) then
              do m2=-spk(ist1,is),spk(ist1,is)-1
                ic=ic+1
                call wavefcr(lradstp,is,ia,ist1,m2,nrcmtmax,wfcr2)
                zfmt(:,1:nrc)=conjg(wfcr1(:,1:nrc,1))*wfcr2(:,1:nrc,1) &
                             +conjg(wfcr1(:,1:nrc,2))*wfcr2(:,1:nrc,2)
                call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr, &
                 zfmt,lmmaxvr,zzero,zrhomt2(:,:,ic),lmmaxvr)
              end do
            end if
          end do
          do ist2=1,nstsv
            if (evalsv(ist2,ikp).gt.efermi) then
! calculate the Coulomb potential
              call zpotclmt(lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
               zrhomt1(:,:,ias,ist2),zvclmt)
!-------------------------------------------!
!     valence-core-valence contribution     !
!-------------------------------------------!
              do ist1=1,nstsv
                if (evalsv(ist1,ikp).lt.efermi) then
                  zt1=zfmtinp(.true.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
                   zrhomt1(:,:,ias,ist1),zvclmt)
                  vnlvv(ist1,ist2)=vnlvv(ist1,ist2)-zt1
                end if
              end do
!----------------------------------------!
!     core-core-valence contribution     !
!----------------------------------------!
              ic=0
              do ist1=1,spnst(is)
                if (spcore(ist1,is)) then
                  do m2=-spk(ist1,is),spk(ist1,is)-1
                    ic=ic+1
                    zt1=zfmtinp(.true.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
                     zrhomt2(:,:,ic),zvclmt)
                    vnlcv(ic,ias,ist2)=vnlcv(ic,ias,ist2)-zt1
                  end do
! end loop over ist1
                end if
              end do
! end loop over ist2
            end if
          end do
! end loops over ist3 and m1
        end do
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr)
deallocate(vgqc,tpgqc,gqc,jlgqr,jlgq0r)
deallocate(sfacgknr,apwalm,evecfv,evecsv,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2,wfcr1,wfcr2)
deallocate(zrhomt1,zrhomt2,zrhoir1)
deallocate(zvclmt,zvclir,zfmt)
return
end subroutine
!EOC

