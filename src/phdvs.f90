
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdvs(p,vsmt1,vsir1)
use modmain
use modphonon
use modstore
implicit none
! arguments
integer, intent(in) :: p
real(8), intent(in) :: vsmt1(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: vsir1(ngtot)
! local variables
integer is,ia,ja,ias,jas
integer ir,irc,i1,i2,i3,i
real(8) v1(3),v2(3),v3(3),t1,t2
complex(8) z1,z2
! automatic arrays
real(8) rflm(lmmaxvr)
complex(8) zflm(lmmaxvr)
! external functions
real(8) rfirvec
external rfirvec
! prefactor
z1=1.d0/(dble(nscph)*deltaph)
! multiply by i for sin-like displacement
if (p.eq.1) z1=z1*zi
!------------------------------!
!     muffin-tin potential     !
!------------------------------!
ias=0
jas=0
do is=1,nspecies
  ja=0
  do ia=1,natoms0(is)
    ias=ias+1
    do i=1,nscph
      ja=ja+1
      jas=jas+1
! important: the muffin-tin potential should have an *explicit* phase exp(iq.r)
      t1=-dot_product(vqc(:,iqph),atposc(:,ja,is))
      z2=z1*cmplx(cos(t1),sin(t1),8)
! loop over radial points
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
! compute the difference between the perturbed and unperturbed potentials
        rflm(:)=vsmt(:,ir,jas)-vsmt1(:,ir,jas)
! convert real potential to a complex spherical harmonic expansion
        call rtozflm(lmaxvr,rflm,zflm)
! add to total
        dvsmt(:,irc,ias)=dvsmt(:,irc,ias)+z2*zflm(:)
! end loop over radial points
      end do
    end do
! end loop over atoms and species
  end do
end do
!--------------------------------!
!     interstitial potential     !
!--------------------------------!
ir=0
do i3=0,ngridg0(3)-1
  v1(3)=dble(i3)/dble(ngridg0(3))
  do i2=0,ngridg0(2)-1
    v1(2)=dble(i2)/dble(ngridg0(2))
    do i1=0,ngridg0(1)-1
      v1(1)=dble(i1)/dble(ngridg0(1))
      ir=ir+1
      call r3mv(avec0,v1,v2)
      do i=1,nscph
        v3(:)=v2(:)+vscph(:,i)
!***** remove vscph
        t1=-dot_product(vqc(:,iqph),v3(:))
        z2=z1*cmplx(cos(t1),sin(t1),8)
        t1=rfirvec(ngridg,ainv,v3,vsir)
        t2=rfirvec(ngridg,ainv,v3,vsir1)
        dvsir(ir)=dvsir(ir)+z2*(t1-t2)
      end do
    end do
  end do
end do
!*********** fix me!
!***** FFT first
return
end subroutine

