
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genwfhk(ik,vhkc,wfmt,wfir,wfhk)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: vhkc(3,nhvec)
complex(8), intent(inout) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfir(ngkmax,nspinor,nstsv)
complex(8), intent(out) :: wfhk(nhvec,nspinor,nstsv)
! local variables
integer igk,ih,ist,ispn
integer is,ia,ias,nrc,irc
integer l,m,lm,i
real(8) v(3),x,hk,tp(2)
real(8) t0,t1,t2
complex(8) zsum1,zsum2
complex(8) zt1,zt2,zt3
! automatic arrays
real(8) jl(0:lmaxvr,nrcmtmax)
real(8) fr1(nrcmtmax),fr2(nrcmtmax),gr(nrcmtmax)
complex(8) ylm(lmmaxvr)
! zero the coefficients
wfhk(:,:,:)=0.d0
!---------------------------!
!     interstitial part     !
!---------------------------!
i=1
do ih=1,nhvec
  v(:)=vhkc(:,ih)
  do igk=i,ngk(1,ik)
    t1=abs(v(1)-vgkc(1,igk,1,ik)) &
      +abs(v(2)-vgkc(2,igk,1,ik)) &
      +abs(v(3)-vgkc(3,igk,1,ik))
    if (t1.lt.epslat) then
      do ist=1,nstsv
        do ispn=1,nspinor
          wfhk(ih,ispn,ist)=wfir(igk,ispn,ist)
        end do
      end do
      if (igk.eq.i) i=igk
    end if
  end do
end do
!-------------------------!
!     muffin-tin part     !
!-------------------------!
t0=fourpi/sqrt(omega)
! remove continuation of interstitial function into muffin-tin
do igk=1,ngk(1,ik)
! generate the spherical harmonics Y_lm(G+k)
  call genylm(lmaxvr,tpgkc(:,igk,1,ik),ylm)
! loop over species
  do is=1,nspecies
    nrc=nrcmt(is)
! generate spherical Bessel functions
    do irc=1,nrc
      x=gkc(igk,1,ik)*rcmt(irc,is)
      call sbessel(lmaxvr,x,jl(:,irc))
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ist=1,nstsv
        do ispn=1,nspinor
          zt1=t0*wfir(igk,ispn,ist)*sfacgk(igk,ias,1,ik)
          lm=0
          do l=0,lmaxvr
            zt2=zt1*zil(l)
            do m=-l,l
              lm=lm+1
              zt3=zt2*conjg(ylm(lm))
              wfmt(lm,1:nrc,ias,ispn,ist)=wfmt(lm,1:nrc,ias,ispn,ist) &
               -zt3*jl(l,1:nrc)
            end do
          end do
        end do
      end do
    end do
  end do
end do
! loop over H+k-vectors
do ih=1,nhvec
! spherical coordinates of H+k
  call sphcrd(vhkc(:,ih),hk,tp)
! generate the spherical harmonics Y_lm(H+k)
  call genylm(lmaxvr,tp,ylm)
  do is=1,nspecies
    nrc=nrcmt(is)
! generate spherical Bessel functions
    do irc=1,nrc
      x=hk*rcmt(irc,is)
      call sbessel(lmaxvr,x,jl(:,irc))
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! conjugate structure factor
      t1=-dot_product(vhkc(:,ih),atposc(:,ia,is))
      zt1=cmplx(cos(t1),sin(t1),8)
      do ist=1,nstsv
        do ispn=1,nspinor
          do irc=1,nrc
            zsum1=0.d0
            lm=0
            do l=0,lmaxvr
              zsum2=0.d0
              do m=-l,l
                lm=lm+1
                zsum2=zsum2+wfmt(lm,irc,ias,ispn,ist)*ylm(lm)
              end do
              zsum1=zsum1+jl(l,irc)*conjg(zil(l))*zsum2
            end do
            zsum1=zsum1*rcmt(irc,is)**2
            fr1(irc)=dble(zsum1)
            fr2(irc)=aimag(zsum1)
          end do
          call fderiv(-1,nrc,rcmt(:,is),fr1,gr)
          t1=gr(nrc)
          call fderiv(-1,nrc,rcmt(:,is),fr2,gr)
          t2=gr(nrc)
          wfhk(ih,ispn,ist)=wfhk(ih,ispn,ist)+t0*zt1*cmplx(t1,t2,8)
        end do
      end do
    end do
  end do
end do
return
end subroutine

