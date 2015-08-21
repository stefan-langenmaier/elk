
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvbmatk(vmt,vir,bmt,bir,wfmt,wfir,vbmat)
!***** remember to multiply by cfunir!!!
use modmain
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(lmmaxvr,nrcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
complex(8), intent(in) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfir(ngtot,nspinor,nstsv)
complex(8), intent(out) :: vbmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,is,ias
integer nrc,nrci,ir,irc
integer lmmax,itp
real(8) t0,t1,t2
complex(8) z1
! allocatable arrays
complex(8), allocatable :: zfmt(:,:,:),zfir(:,:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
! zero the matrix elements
vbmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(zfmt(lmmaxvr,nrcmtmax,nspinor))
do jst=1,nstsv
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmtinr(is)
! apply local potential plus magnetic field to spinor wavefunction
    lmmax=lmmaxinr
    do irc=1,nrc
      do itp=1,lmmax
        t1=vmt(itp,irc,ias)
        t2=bmt(itp,irc,ias,ndmag)
        zfmt(itp,irc,1)=(t1+t2)*wfmt(itp,irc,ias,1,jst)
        zfmt(itp,irc,2)=(t1-t2)*wfmt(itp,irc,ias,2,jst)
      end do
      if (ncmag) then
        do itp=1,lmmax
          z1=cmplx(bmt(itp,irc,ias,1),bmt(itp,irc,ias,2),8)
          zfmt(itp,irc,1)=zfmt(itp,irc,1)+conjg(z1)*wfmt(itp,irc,ias,2,jst)
          zfmt(itp,irc,2)=zfmt(itp,irc,2)+z1*wfmt(itp,irc,ias,1,jst)
        end do
      end if
      if (irc.eq.nrci) lmmax=lmmaxvr
    end do
    do ist=1,jst
! compute inner product (functions are in spherical coordinates)
      do ispn=1,nspinor
        z1=zfmtinp(.false.,nrc,nrci,rcmt(:,is),wfmt(:,:,ias,ispn,ist), &
         zfmt(:,:,ispn))
        vbmat(ist,jst)=vbmat(ist,jst)+z1
      end do
    end do
  end do
end do
deallocate(zfmt)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(zfir(ngtot,nspinor))
t0=omega/dble(ngtot)
do jst=1,nstsv
! apply local potential and magnetic field to spinor wavefuntion
  do ir=1,ngtot
    t1=vir(ir)
    t2=bir(ir,ndmag)
    zfir(ir,1)=(t1+t2)*wfir(ir,1,jst)
    zfir(ir,2)=(t1-t2)*wfir(ir,2,jst)
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
      vbmat(ist,jst)=vbmat(ist,jst)+t0*z1
    end do
  end do
end do
deallocate(zfir)
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vbmat(ist,jst)=conjg(vbmat(jst,ist))
  end do
end do
return
end subroutine
