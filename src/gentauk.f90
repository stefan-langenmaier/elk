
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentauk(ik,taumt,tauir)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(inout) :: taumt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(inout) :: tauir(ngrtot)
! local variables
integer ist,is,ias,ispn,jspn
integer ir,irc,itp,igk,ifg
real(8) wo,t0,t1
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: zfir(:)
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngkmax,nspinor,nstsv))
allocate(zfir(ngrtot))
! generate the second-variational wavefunctions for all states
call genwfsvp(.false.,.true.,vkl(:,ik),wfmt,ngkmax,wfir)
! loop over states
do ist=1,nstsv
  wo=wkpt(ik)*occsv(ist,ik)
  if (abs(wo).lt.epsocc) cycle
  t0=wo*evalsv(ist,ik)
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
  do ias=1,natmtot
    is=idxis(ias)
    do ispn=1,nspinor
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
        do itp=1,lmmaxvr
          zt1=wfmt(itp,irc,ias,ispn,ist)
          t1=dble(zt1)**2+aimag(zt1)**2
          taumt(itp,ir,ias)=taumt(itp,ir,ias)+t0*t1
        end do
      end do
    end do
  end do
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
  t0=t0/omega
  do ispn=1,nspinor
    if (spinsprl) then
      jspn=ispn
    else
      jspn=1
    end if
    zfir(:)=0.d0
    do igk=1,ngk(jspn,ik)
      ifg=igfft(igkig(igk,jspn,ik))
      zfir(ifg)=wfir(igk,ispn,ist)
    end do
! Fourier transform to real-space
    call zfftifc(3,ngrid,1,zfir)
    do ir=1,ngrtot
      t1=dble(zfir(ir))**2+aimag(zfir(ir))**2
      tauir(ir)=tauir(ir)+t0*t1
    end do
  end do
! end loop over states
end do
deallocate(wfmt,wfir,zfir)
return
end subroutine

