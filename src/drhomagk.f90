
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
integer nst,ist,jst,is,ias
integer nr,nrci,ir,irc
integer lmmax,itp
real(8) t0
complex(8) z1,z2,z3,z4,z5,z6
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: dwfmt(:,:,:,:,:),dwfir(:,:,:)
! count and index the occupied states
nst=0
do ist=1,nstsv
  if (abs(occsvp(ist)).gt.epsocc) then
    nst=nst+1
    idx(nst)=ist
  end if
end do
! generate the wavefunctions
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nst))
allocate(wfir(ngtot,nspinor,nst))
call genwfsv(.false.,.false.,nst,idx,ngp,igpig,apwalm,evecfv,evecsv,wfmt, &
 ngtot,wfir)
! generate the wavefunction derivatives
allocate(dwfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nst))
allocate(dwfir(ngtot,nspinor,nst))
call gendwfsv(.false.,.false.,nst,idx,ngp,ngpq,igpqig,apwalmq,dapwalm,evecfv, &
 devecfv,evecsv,devecsv,dwfmt,ngtot,dwfir)
! loop over occupied states
do ist=1,nst
  jst=idx(ist)
  t0=2.d0*wkptnr*occsvp(jst)
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
        lmmax=lmmaxinr
        irc=0
        do ir=1,nr,lradstp
          irc=irc+1
          do itp=1,lmmax
            z1=conjg(wfmt(itp,irc,ias,1,jst))
            z2=conjg(wfmt(itp,irc,ias,2,jst))
            z3=dwfmt(itp,irc,ias,1,jst)
            z4=dwfmt(itp,irc,ias,2,jst)
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
          if (irc.eq.nrci) lmmax=lmmaxvr
        end do
      else
! collinear
        lmmax=lmmaxinr
        irc=0
        do ir=1,nr,lradstp
          irc=irc+1
          do itp=1,lmmax
            z1=conjg(wfmt(itp,irc,ias,1,jst))*dwfmt(itp,irc,ias,1,jst)
            z2=conjg(wfmt(itp,irc,ias,2,jst))*dwfmt(itp,irc,ias,2,jst)
            drhomt(itp,ir,ias)=drhomt(itp,ir,ias)+t0*(z1+z2)
            dmagmt(itp,ir,ias,1)=dmagmt(itp,ir,ias,1)+t0*(z1-z2)
          end do
          if (irc.eq.nrci) lmmax=lmmaxvr
        end do
      end if
    else
! spin-unpolarised
      lmmax=lmmaxinr
      irc=0
      do ir=1,nr,lradstp
        irc=irc+1
        drhomt(1:lmmax,ir,ias)=drhomt(1:lmmax,ir,ias) &
         +t0*conjg(wfmt(1:lmmax,irc,ias,1,jst))*dwfmt(1:lmmax,irc,ias,1,jst)
        if (irc.eq.nrci) lmmax=lmmaxvr
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
        z1=conjg(wfir(ir,1,jst))
        z2=conjg(wfir(ir,2,jst))
        z3=dwfir(ir,1,jst)
        z4=dwfir(ir,2,jst)
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
        z1=conjg(wfir(ir,1,jst))*dwfir(ir,1,jst)
        z2=conjg(wfir(ir,2,jst))*dwfir(ir,2,jst)
        drhoir(ir)=drhoir(ir)+t0*(z1+z2)
        dmagir(ir,1)=dmagir(ir,1)+t0*(z1-z2)
      end do
    end if
  else
! spin-unpolarised
    do ir=1,ngtot
      drhoir(ir)=drhoir(ir)+t0*conjg(wfir(ir,1,jst))*dwfir(ir,1,jst)
    end do
  end if
!$OMP END CRITICAL
! end loop over states
end do
deallocate(wfmt,wfir,dwfmt,dwfir)
return
end subroutine

