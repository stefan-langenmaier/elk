
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdwfsv(ik,dwfmt,dwfir)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(out) :: dwfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: dwfir(ngrtot,nspinor,nstsv)
! local variables
integer jk,jkq,ist,jst
integer igk,igkq,iv(3),ifg
integer ispn0,ispn1
integer ispn,jspn,is,ia,ias
integer nrc,irc,ir,l,m,lm
real(8) eij,x,t0
complex(8) zt1,zt2,zt3
! automatic arrays
complex(8) ylm(lmmaxvr)
! allocatable arrays
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: wfpw(:,:,:),wfpwq(:,:,:)
complex(8), allocatable :: wfpwh(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:),zv(:,:)
complex(8), allocatable :: gzfmt(:,:,:)
! external functions
complex(8) zdotc
external zdotc
allocate(jl(0:lmaxvr,nrcmtmax))
allocate(wfpw(ngkmax,nspinor,nstsv),wfpwq(ngkmax,nspinor,nstsv))
allocate(wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor),zv(ngkmax,nspinor))
allocate(gzfmt(lmmaxvr,nrcmtmax,3))
! equivalent reduced k-point
jk=ikmap(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! equivalent reduced k+q-point
iv(:)=ivk(:,ik)+ivq(:,iqph)
iv(:)=modulo(iv(:),ngridk(:))
jkq=ikmap(iv(1),iv(2),iv(3))
! generate the low and high plane wave wavefunctions at k
call genwfpw(.true.,vkl(:,ik),ngk(1,ik),igkig(:,1,ik),vgkl(:,:,1,ik), &
 gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),wfpw,wfpwh)
! generate the low plane wave wavefunction at k+q
call genwfpw(.false.,vkql(:,ik),ngkq(1,ik),igkqig(:,1,ik),vgkql(:,:,1,ik), &
 gkqc(:,1,ik),tpgkqc(:,:,1,ik),sfacgkq(:,:,1,ik),wfpwq,wfpwh)
! loop over states
do ist=1,nstsv
! Fourier transform wavefunction to real-space
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    wfir(:,ispn)=0.d0
    do igk=1,ngk(jspn,ik)
      ifg=igfft(igkig(igk,jspn,ik))
      wfir(ifg,ispn)=wfpw(igk,ispn,ist)
    end do
    call zfftifc(3,ngrid,1,wfir(:,ispn))
  end do
! multiply wavefunction by effective potential derivative
  if (spinpol) then
! spin polarised
    if (ncmag) then
! non-collinear
      do ir=1,ngrtot
        zt1=wfir(ir,1); zt2=wfir(ir,2)
        wfir(ir,1)=dveffpw(1,1,ir)*zt1+dveffpw(1,2,ir)*zt2
        wfir(ir,2)=dveffpw(2,1,ir)*zt1+dveffpw(2,2,ir)*zt2
      end do
    else
! collinear
      do ir=1,ngrtot
        wfir(ir,1)=dveffpw(1,1,ir)*wfir(ir,1)
        wfir(ir,2)=dveffpw(2,2,ir)*wfir(ir,2)
      end do
    end if
  else
! spin-unpolarised
    do ir=1,ngrtot
      wfir(ir,1)=dveffpw(1,1,ir)*wfir(ir,1)
    end do
  end if
! Fourier transform back to G+k+q-space
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    call zfftifc(3,ngrid,-1,wfir(:,ispn))
    do igkq=1,ngkq(jspn,ik)
      ifg=igfft(igkqig(igkq,jspn,ik))
      zv(igkq,ispn)=wfir(ifg,ispn)
    end do
  end do
  dwfir(:,:,ist)=0.d0
  do jst=1,nstsv
    eij=evalsv(ist,jk)-evalsv(jst,jkq)
    if (abs(eij).gt.1.d-8) then
      zt1=0.d0
      do ispn=1,nspinor
        jspn=jspnfv(ispn)
        zt1=zt1+zdotc(ngkq(jspn,ik),wfpwq(:,ispn,jst),1,zv(:,ispn),1)
      end do
      zt1=zt1/eij
      do ispn=1,nspinor
        jspn=jspnfv(ispn)
        call zaxpy(ngkq(jspn,ik),zt1,wfpwq(:,ispn,jst),1,dwfir(:,ispn,ist),1)
      end do
    end if
  end do
end do
! plane wave extension into muffin-tins
t0=fourpi/sqrt(omega)
dwfmt(:,:,:,:,:)=0.d0
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
  do igkq=1,ngkq(jspn,ik)
! generate the spherical harmonics Y_lm(G+k+q)
    call genylm(lmaxvr,tpgkqc(:,igkq,jspn,ik),ylm)
! loop over species
    do is=1,nspecies
      nrc=nrcmt(is)
! generate spherical Bessel functions
      do irc=1,nrc
        x=gkqc(igkq,jspn,ik)*rcmt(irc,is)
        call sbessel(lmaxvr,x,jl(:,irc))
      end do
! loop over atoms
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        zt1=t0*sfacgkq(igkq,ias,jspn,ik)
! loop over states
        do ist=1,nstsv
          do ispn=ispn0,ispn1
            zt2=zt1*dwfir(igkq,ispn,ist)
            do irc=1,nrc
              lm=0
              do l=0,lmaxvr
                zt3=jl(l,irc)*zil(l)*zt2
                do m=-l,l
                  lm=lm+1
                  dwfmt(lm,irc,ias,ispn,ist)=dwfmt(lm,irc,ias,ispn,ist) &
                   +zt3*conjg(ylm(lm))
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do
! Fourier transform interstitial part to real-space
t0=1.d0/sqrt(omega)
do ist=1,nstsv
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
    do igkq=1,ngkq(jspn,ik)
      zv(igkq,ispn)=t0*dwfir(igkq,ispn,ist)
    end do
    dwfir(:,ispn,ist)=0.d0
    do igkq=1,ngkq(jspn,ik)
      ifg=igfft(igkqig(igkq,jspn,ik))
      dwfir(ifg,ispn,ist)=zv(igkq,ispn)
    end do
    call zfftifc(3,ngrid,1,dwfir(:,ispn,ist))
  end do
end do
deallocate(jl,wfpw,wfpwq,wfir,zv)
! subtract the gradient of the high plane wave wavefunction
do ist=1,nstsv
  do ispn=1,nspinor
    ias=idxas(iaph,isph)
    nrc=nrcmt(isph)
    call gradzfmt(lmaxvr,nrc,rcmt(:,isph),lmmaxvr,nrcmtmax, &
     wfpwh(:,:,ias,ispn,ist),gzfmt)
    do irc=1,nrc
      dwfmt(:,irc,ias,ispn,ist)=dwfmt(:,irc,ias,ispn,ist)-gzfmt(:,irc,ipph)
    end do
  end do
end do
deallocate(wfpwh,gzfmt)
return
end subroutine

