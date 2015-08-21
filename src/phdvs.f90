
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
integer nrc,nrci,ir,irc
integer i1,i2,i3,i
real(8) v1(3),v2(3),v3(3),t1,t2
complex(8) z1,z2
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
complex(8), allocatable :: zfmt(:,:)
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
allocate(rfmt(lmmaxvr,nrcmtmax),zfmt(lmmaxvr,nrcmtmax))
ias=0
jas=0
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
  ja=0
  do ia=1,natoms0(is)
    ias=ias+1
    do i=1,nscph
      ja=ja+1
      jas=jas+1
! compute the difference between the perturbed and unperturbed potentials
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
        rfmt(:,irc)=vsmt(:,ir,jas)-vsmt1(:,ir,jas)
      end do
! convert real potential difference to a complex spherical harmonic expansion
      call rtozfmt(nrc,nrci,1,rfmt,1,zfmt)
! important: the muffin-tin potential should have an *explicit* phase exp(iq.r)
      t1=-dot_product(vqc(:,iqph),atposc(:,ja,is))
      z2=z1*cmplx(cos(t1),sin(t1),8)
! add to total
      call zfmtadd(nrc,nrci,z2,zfmt,dvsmt(:,:,ias))
    end do
! end loop over atoms and species
  end do
end do
deallocate(rfmt,zfmt)
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

