
! Copyright (C) 2010 Alexey I. Baranov.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zftrf
! !INTERFACE:
subroutine zftrf(rfmt,rfir,zfh)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   rfmt : real muffin-tin function (in,real(lmmaxvr,nrmtmax,natmtot))
!   rfir : real interstitial function (in,real(ngrtot))
!   zfh  : Fourier expansion coefficients of the real-space function
!          (out,complex(nhv))
! !DESCRIPTION:
!   Given a real function, $f({\bf r})$, this routine calculates its complex
!   Fourier expansion coefficients:
!   $$ f({\bf H})=\frac{1}{\Omega}\int d^3r\,f({\bf r})\tilde{\Theta}({\bf r})
!    e^{-i{\bf H}\cdot{\bf r}}
!    +\frac{4\pi}{\Omega}\sum_{\alpha}e^{-i{\bf H}\cdot{\bf R}_{\alpha}}
!    \sum_{lm}(-i)^l Y_{lm}(\hat{\bf H})
!    \int_{0}^{R_{\alpha}}dr\,r^2 j_{l}(|{\bf H}|r)f_{lm}^{\alpha}(r), $$
!   where $\tilde{\Theta}$ is the smooth characteristic function of the
!   interstitial region, $\Omega$ is the unit cell volume and $R_{\alpha}$ is
!   the muffin-tin radius of atom $\alpha$. The $\bf H$-vectors are stored in
!   the global array {\tt ivh}. See also {\tt genhvec}.
!
! !REVISION HISTORY:
!   Created July 2010 (Alexey I. Baranov)
!   Modified, November 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: rfir(ngrtot)
complex(8), intent(out) :: zfh(nhvec)
! local variables
integer is,ia,ias
integer nrc,irc,ir
integer ih,ig,l,m,lm
real(8) vhl(3),vhc(3),hc
real(8) tp(2),x,t1,t2
complex(8) zsum1,zsum2
complex(8) zt1,zt2,zt3
! automatic arrays
real(8) jl(0:lmaxvr,nrcmtmax)
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
real(8) gr(nrcmtmax),cf(4,nrcmtmax)
complex(8) ylm(lmmaxvr)
! allocatable arrays
complex(8), allocatable :: zfft(:)
complex(8), allocatable :: zfmt(:,:,:)
allocate(zfft(ngrtot))
allocate(zfmt(lmmaxvr,nrcmtmax,natmtot))
! zero the coefficients
zfh(:)=0.d0
!---------------------------!
!     interstitial part     !
!---------------------------!
! Fourier transform to G-space
zfft(:)=rfir(:)
call zfftifc(3,ngrid,-1,zfft)
! find coefficients for all required input vectors
do ih=1,nhvec
  if ((ivh(1,ih).ge.intgv(1,1)).and.(ivh(1,ih).le.intgv(1,2)).and. &
      (ivh(2,ih).ge.intgv(2,1)).and.(ivh(2,ih).le.intgv(2,2)).and. &
      (ivh(3,ih).ge.intgv(3,1)).and.(ivh(3,ih).le.intgv(3,2))) then
    ig=ivgig(ivh(1,ih),ivh(2,ih),ivh(3,ih))
    zfh(ih)=zfft(igfft(ig))
  end if
end do
!-------------------------!
!     muffin-tin part     !
!-------------------------!
! convert function from real to complex spherical harmonic expansion
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      call rtozflm(lmaxvr,rfmt(:,ir,ias),zfmt(:,irc,ias))
    end do
  end do
end do
! remove continuation of interstitial function into muffin-tin
do ig=1,ngvec
  do is=1,nspecies
    nrc=nrcmt(is)
! generate spherical Bessel functions
    do irc=1,nrc
      x=gc(ig)*rcmt(irc,is)
      call sbessel(lmaxvr,x,jl(:,irc))
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      zt1=zfft(igfft(ig))*fourpi*sfacg(ig,ias)
      lm=0
      do l=0,lmaxvr
        zt2=zt1*zil(l)
        do m=-l,l
          lm=lm+1
          zt3=zt2*conjg(ylmg(lm,ig))
          zfmt(lm,1:nrc,ias)=zfmt(lm,1:nrc,ias)-zt3*jl(l,1:nrc)
        end do
      end do
    end do
  end do
end do
! loop over input H-vectors
do ih=1,nhvec
  vhl(:)=dble(ivh(:,ih))
  call r3mv(bvec,vhl,vhc)
  call sphcrd(vhc,hc,tp)
  call genylm(lmaxvr,tp,ylm)
  do is=1,nspecies
    nrc=nrcmt(is)
! generate spherical Bessel functions
    do irc=1,nrc
      x=hc*rcmt(irc,is)
      call sbessel(lmaxvr,x,jl(:,irc))
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! structure factor
      t1=dot_product(vhc,atposc(:,ia,is))
      zt1=cmplx(cos(t1),sin(t1),8)
      do irc=1,nrc
        zsum1=0.d0
        lm=0
        do l=0,lmaxvr
          zsum2=0.d0
          do m=-l,l
            lm=lm+1
            zsum2=zsum2+zfmt(lm,irc,ias)*ylm(lm)
          end do
          zsum1=zsum1+jl(l,irc)*conjg(zil(l))*zsum2
        end do
        zsum1=zsum1*rcmt(irc,is)**2
        fr1(irc)=dble(zsum1)
        fr2(irc)=aimag(zsum1)
      end do
      call fderiv(-1,nrc,rcmt(:,is),fr1,gr,cf)
      t1=gr(nrc)
      call fderiv(-1,nrc,rcmt(:,is),fr2,gr,cf)
      t2=gr(nrc)
      zfh(ih)=zfh(ih)+(fourpi/omega)*conjg(zt1)*cmplx(t1,t2,8)
    end do
  end do
end do
deallocate(zfft,zfmt)
return
end subroutine
! EOC

