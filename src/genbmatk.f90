
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genbmatk(bmt,bir,wfmt,wfir,bmat)
! calculates the magnetic field matrix elements
use modmain
implicit none
! arguments
real(8), intent(in) :: bmt(lmmaxvr,nrcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfir(ngtot,nspinor,nstsv)
complex(8), intent(out) :: bmat(nstsv,nstsv)
! local variables
integer ispn,ist,jst
integer is,ias,nrc,ir,irc
real(8) t1
complex(8) z1
! automatic arrays
complex(8) zflm(lmmaxvr)
! allocatable arrays
complex(8), allocatable :: zfmt(:,:,:),zfir(:,:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
allocate(zfmt(lmmaxvr,nrcmtmax,nspinor))
allocate(zfir(ngtot,nspinor))
! zero the matrix elements
bmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do jst=1,nstsv
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
! apply magnetic field to spinor wavefunction
    do irc=1,nrc
      zfmt(:,irc,1)=bmt(:,irc,ias,ndmag)*wfmt(:,irc,ias,1,jst)
      zfmt(:,irc,2)=-bmt(:,irc,ias,ndmag)*wfmt(:,irc,ias,2,jst)
      if (ncmag) then
        zflm(:)=cmplx(bmt(:,irc,ias,1),bmt(:,irc,ias,2),8)
        zfmt(:,irc,1)=zfmt(:,irc,1)+conjg(zflm(:))*wfmt(:,irc,ias,2,jst)
        zfmt(:,irc,2)=zfmt(:,irc,2)+zflm(:)*wfmt(:,irc,ias,1,jst)
      end if
    end do
    do ist=1,jst
! compute inner product (functions are in spherical coordinates)
      do ispn=1,nspinor
        z1=zfmtinp(.false.,lmmaxvr,nrc,rcmt(:,is),lmmaxvr, &
         wfmt(:,:,ias,ispn,ist),zfmt(:,:,ispn))
        bmat(ist,jst)=bmat(ist,jst)+z1
      end do
    end do
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
t1=omega/dble(ngtot)
do jst=1,nstsv
! apply magnetic field to spinor wavefunction
  do ir=1,ngtot
    zfir(ir,1)=bir(ir,ndmag)*wfir(ir,1,jst)
    zfir(ir,2)=-bir(ir,ndmag)*wfir(ir,2,jst)
  end do
  if (ncmag) then
    do ir=1,ngtot
      z1=cmplx(bir(ir,1),bir(ir,2),8)
      zfir(ir,1)=zfir(ir,1)+conjg(z1)*wfir(ir,2,jst)
      zfir(ir,2)=zfir(ir,2)+z1*wfir(ir,1,jst)
    end do
  end if
  do ist=1,jst
    do ispn=1,nspinor
      z1=zdotc(ngtot,wfir(:,ispn,ist),1,zfir(:,ispn),1)
      bmat(ist,jst)=bmat(ist,jst)+t1*z1
    end do
  end do
end do
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    bmat(ist,jst)=conjg(bmat(jst,ist))
  end do
end do
deallocate(zfmt,zfir)
return
end subroutine

