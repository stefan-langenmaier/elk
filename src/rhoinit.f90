
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhoinit
! !INTERFACE:
subroutine rhoinit
! !USES:
use modmain
! !DESCRIPTION:
!   Initialises the crystal charge density. Inside the muffin-tins it is set to
!   the spherical atomic density. In the interstitial region it is taken to be
!   constant such that the total charge is correct. Requires that the atomic
!   densities have already been calculated.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: lmax=1
integer lmmax,l,m,lm,ir,irc
integer is,ia,ias,ig,ifg
real(8) x,t1,t2
complex(8) zt1,zt2,zt3
! automatic arrays
real(8) fr(spnrmax),gr(spnrmax),cf(4,spnrmax)
! allocatable arrays
real(8), allocatable :: jlgr(:,:)
real(8), allocatable :: th(:)
real(8), allocatable :: ffacg(:)
complex(8), allocatable :: zfmt(:,:)
complex(8), allocatable :: zfft(:)
! external functions
real(8) erf
external erf
lmmax=(lmax+1)**2
! allocate local arrays
allocate(jlgr(0:lmax,nrcmtmax))
allocate(th(spnrmax))
allocate(ffacg(ngvec))
allocate(zfmt(lmmax,nrcmtmax))
allocate(zfft(ngrtot))
! zero the charge density and magnetisation arrays
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  magir(:,:)=0.d0
end if
! compute the superposition of all the atomic density tails
zfft(:)=0.d0
do is=1,nspecies
! generate smooth step function
  do ir=1,spnr(is)
    x=0.5d0*gmaxvr*(spr(ir,is)-rmt(is))
    th(ir)=0.5d0*(1.d0+erf(x))
  end do
  do ig=1,ngvec
    do ir=1,spnr(is)
! spherical bessel function j_0(x)
      x=gc(ig)*spr(ir,is)
      if (x.gt.1.d-6) then
        t1=sin(x)/x
      else
        t1=1.d0
      end if
      fr(ir)=t1*th(ir)*sprho(ir,is)*spr(ir,is)**2
    end do
    call fderiv(-1,spnr(is),spr(:,is),fr,gr,cf)
    ffacg(ig)=(fourpi/omega)*gr(spnr(is))
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ig=1,ngvec
      ifg=igfft(ig)
      zfft(ifg)=zfft(ifg)+ffacg(ig)*conjg(sfacg(ig,ias))
    end do
  end do
end do
! compute the tails in each muffin-tin
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zfmt(:,:)=0.d0
    do ig=1,ngvec
      ifg=igfft(ig)
      do irc=1,nrcmt(is)
        x=gc(ig)*rcmt(irc,is)
        call sbessel(lmax,x,jlgr(:,irc))
      end do
      zt1=fourpi*zfft(ifg)*sfacg(ig,ias)
      lm=0
      do l=0,lmax
        zt2=zt1*zil(l)
        do m=-l,l
          lm=lm+1
          zt3=zt2*conjg(ylmg(lm,ig))
          do irc=1,nrcmt(is)
            zfmt(lm,irc)=zfmt(lm,irc)+jlgr(l,irc)*zt3
          end do
        end do
      end do
    end do
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      call ztorflm(lmax,zfmt(:,irc),rhomt(:,ir,ias))
    end do
  end do
end do
! convert the density from a coarse to a fine radial mesh
call rfmtctof(rhomt)
! add the atomic charge density and the excess charge in each muffin-tin
t1=chgexs/omega
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      t2=(t1+sprho(ir,is))/y00
      rhomt(1,ir,ias)=rhomt(1,ir,ias)+t2
    end do
  end do
end do
! interstitial density determined from the atomic tails and excess charge
call zfftifc(3,ngrid,1,zfft)
do ir=1,ngrtot
  rhoir(ir)=dble(zfft(ir))+t1
end do
! compute the total charge
call charge
! normalise the density
call rhonorm
deallocate(jlgr,th,ffacg,zfmt,zfft)
return
end subroutine
!EOC

