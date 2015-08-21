
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine drhomagk(ngp,ngpq,igpig,igpqig,occsvp,apwalm,apwalmq,dapwalm, &
 evecfv,devecfv,evecsv,devecsv)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv),ngpq(nspnfv)
integer, intent(in) :: igpig(ngkmax,nspnfv),igpqig(ngkmax,nspnfv)
real(8), intent(in) :: occsvp(nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: apwalmq(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: dapwalm(ngkmax,apwordmax,lmmaxapw,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: devecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv),devecsv(nstsv,nstsv)
! local variables
integer ist,is,ias
integer nr,nrci,ir,irc
integer lmmax,itp
real(8) t0
complex(8) z1,z2,z3,z4,z5,z6
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: dwfmt(:,:,:,:,:),dwfir(:,:,:)
! generate the wavefunctions
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngtot,nspinor,nstsv))
call genwfsv(.false.,.false.,.true.,ngp,igpig,occsvp,apwalm,evecfv,evecsv, &
 wfmt,ngtot,wfir)
! generate the wavefunction derivatives
allocate(dwfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(dwfir(ngtot,nspinor,nstsv))
call gendwfsv(.false.,.false.,.true.,ngp,ngpq,igpqig,occsvp,apwalmq,dapwalm, &
 evecfv,devecfv,evecsv,devecsv,dwfmt,ngtot,dwfir)
! loop over states
do ist=1,nstsv
  if (abs(occsvp(ist)).lt.epsocc) cycle
  t0=2.d0*wkptnr*occsvp(ist)
!----------------------------------------------!
!     muffin-tin density and magnetisation     !
!----------------------------------------------!
  do ias=1,natmtot
    is=idxis(ias)
    nr=nrmt(is)
    nrci=nrcmtinr(is)
!$OMP CRITICAL
    if (spinpol) then
! spin-polarised
      if (ncmag) then
! non-collinear
        irc=0
        do ir=1,nr,lradstp
          irc=irc+1
          if (irc.le.nrci) then
            lmmax=lmmaxinr
          else
            lmmax=lmmaxvr
          end if
          do itp=1,lmmax
            z1=conjg(wfmt(itp,irc,ias,1,ist))
            z2=conjg(wfmt(itp,irc,ias,2,ist))
            z3=dwfmt(itp,irc,ias,1,ist)
            z4=dwfmt(itp,irc,ias,2,ist)
            z5=z1*z3
            z6=z2*z4
            drhomt(itp,ir,ias)=drhomt(itp,ir,ias)+t0*(z5+z6)
            dmagmt(itp,ir,ias,3)=dmagmt(itp,ir,ias,3)+t0*(z5-z6)
            z5=z1*z4
            z6=z2*z3
            dmagmt(itp,ir,ias,1)=dmagmt(itp,ir,ias,1)+t0*(z5+z6)
            z5=z5-z6
            dmagmt(itp,ir,ias,2)=dmagmt(itp,ir,ias,2) &
             +t0*cmplx(aimag(z5),-dble(z5),8)
          end do
        end do
      else
! collinear
        irc=0
        do ir=1,nr,lradstp
          irc=irc+1
          if (irc.le.nrci) then
            lmmax=lmmaxinr
          else
            lmmax=lmmaxvr
          end if
          do itp=1,lmmax
            z1=conjg(wfmt(itp,irc,ias,1,ist))*dwfmt(itp,irc,ias,1,ist)
            z2=conjg(wfmt(itp,irc,ias,2,ist))*dwfmt(itp,irc,ias,2,ist)
            drhomt(itp,ir,ias)=drhomt(itp,ir,ias)+t0*(z1+z2)
            dmagmt(itp,ir,ias,1)=dmagmt(itp,ir,ias,1)+t0*(z1-z2)
          end do
        end do
      end if
    else
! spin-unpolarised
      irc=0
      do ir=1,nr,lradstp
        irc=irc+1
        if (irc.le.nrci) then
          lmmax=lmmaxinr
        else
          lmmax=lmmaxvr
        end if
        drhomt(1:lmmax,ir,ias)=drhomt(1:lmmax,ir,ias) &
         +t0*conjg(wfmt(1:lmmax,irc,ias,1,ist))*dwfmt(1:lmmax,irc,ias,1,ist)
      end do
    end if
!$OMP END CRITICAL
  end do
!------------------------------------------------!
!     interstitial density and magnetisation     !
!------------------------------------------------!
!$OMP CRITICAL
  if (spinpol) then
! spin-polarised
    if (ncmag) then
! non-collinear
      do ir=1,ngtot
        z1=conjg(wfir(ir,1,ist))
        z2=conjg(wfir(ir,2,ist))
        z3=dwfir(ir,1,ist)
        z4=dwfir(ir,2,ist)
        z5=z1*z3
        z6=z2*z4
        drhoir(ir)=drhoir(ir)+t0*(z5+z6)
        dmagir(ir,3)=dmagir(ir,3)+t0*(z5-z6)
        z5=z1*z4
        z6=z2*z3
        dmagir(ir,1)=dmagir(ir,1)+t0*(z5+z6)
        z5=z5-z6
        dmagir(ir,2)=dmagir(ir,2)+t0*cmplx(aimag(z5),-dble(z5),8)
      end do
    else
! collinear
      do ir=1,ngtot
        z1=conjg(wfir(ir,1,ist))*dwfir(ir,1,ist)
        z2=conjg(wfir(ir,2,ist))*dwfir(ir,2,ist)
        drhoir(ir)=drhoir(ir)+t0*(z1+z2)
        dmagir(ir,1)=dmagir(ir,1)+t0*(z1-z2)
      end do
    end if
  else
! spin-unpolarised
    do ir=1,ngtot
      drhoir(ir)=drhoir(ir)+t0*conjg(wfir(ir,1,ist))*dwfir(ir,1,ist)
    end do
  end if
!$OMP END CRITICAL
! end loop over states
end do
deallocate(wfmt,wfir,dwfmt,dwfir)
return
end subroutine

