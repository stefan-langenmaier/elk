
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmatk(vmt,vir,wfmt,wfir,vmat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot)
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfir(ngtot,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,is,ias
integer nrc,nrci,irc
real(8) t1
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zfmt(:,:),zfir(:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(zfmt(lmmaxvr,nrcmtmax))
do jst=1,nstsv
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmtinr(is)
    do ispn=1,nspinor
! apply potential to wavefunction
      do irc=1,nrci
        zfmt(1:lmmaxinr,irc)=vmt(1:lmmaxinr,irc,ias) &
         *wfmt(1:lmmaxinr,irc,ias,ispn,jst)
      end do
      do irc=nrci+1,nrc
        zfmt(:,irc)=vmt(:,irc,ias)*wfmt(:,irc,ias,ispn,jst)
      end do
      do ist=1,jst
! compute inner product (functions are in spherical coordinates)
        z1=zfmtinp(.false.,nrc,nrci,rcmt(:,is),wfmt(:,:,ias,ispn,ist),zfmt)
        vmat(ist,jst)=vmat(ist,jst)+z1
      end do
    end do
  end do
end do
deallocate(zfmt)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(zfir(ngtot))
t1=omega/dble(ngtot)
do jst=1,nstsv
  do ispn=1,nspinor
! apply potential to wavefunction
    zfir(:)=vir(:)*wfir(:,ispn,jst)
    do ist=1,jst
      z1=zdotc(ngtot,wfir(:,ispn,ist),1,zfir,1)
      vmat(ist,jst)=vmat(ist,jst)+t1*z1
    end do
  end do
end do
deallocate(zfir)
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vmat(ist,jst)=conjg(vmat(jst,ist))
  end do
end do
return
end subroutine

