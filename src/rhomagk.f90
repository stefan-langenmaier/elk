
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhomagk
! !INTERFACE:
subroutine rhomagk(ik,evecfv,evecsv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Generates the partial valence charge density from the eigenvectors at
!   $k$-point {\tt ik}. In the muffin-tin region, the wavefunction is obtained
!   in terms of its $(l,m)$-components from both the APW and local-orbital
!   functions. Using a backward spherical harmonic transform (SHT), the
!   wavefunction is converted to real-space and the density obtained from its
!   modulus squared. This density is then accumulated in the global variable
!   {\tt rhomt}. A similar proccess is used for the intersitial density in which
!   the wavefunction in real-space is obtained from a Fourier transform of the
!   sum of APW functions. The interstitial density is added to the global array
!   {\tt rhoir}. See routines {\tt wavefmt}, {\tt genshtmat} and {\tt eveqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Removed conversion to spherical harmonics, January 2009 (JKD)
!   Partially de-phased the muffin-tin magnetisation for spin-spirals,
!    February 2009 (FC, FB & LN)
!   Optimisations, July 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer ispn,jspn,ist,is,ia,ias
integer nr,nrc,nrci,ir,irc
integer lmmax,itp,igk,ifg,i,j
real(8) t0,t1,t2,t3,t4
real(8) ts0,ts1
complex(8) zq(2),z1,z2
! automatic arrays
logical done(nstfv,nspnfv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfmt2(:,:)
complex(8), allocatable :: wfmt3(:,:,:),wfir(:,:)
call timesec(ts0)
!----------------------------------------------!
!     muffin-tin density and magnetisation     !
!----------------------------------------------!
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
if (tevecsv) allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmaxvr,nrcmtmax))
allocate(wfmt3(lmmaxvr,nrcmtmax,nspinor))
! find the matching coefficients
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ispn=1,nspnfv
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
!$OMP END DO
!$OMP END PARALLEL
! loop over atoms
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  nr=nrmt(is)
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
! de-phasing factor for spin-spirals
  if (spinsprl.and.ssdph) then
    t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
    zq(1)=cmplx(cos(t1),sin(t1),8)
    zq(2)=conjg(zq(1))
  end if
  done(:,:)=.false.
  do j=1,nstsv
    if (abs(occsv(j,ik)).lt.epsocc) cycle
    t0=wkpt(ik)*occsv(j,ik)
    t4=2.d0*t0
    if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
      i=0
      do ispn=1,nspinor
        jspn=jspnfv(ispn)
        wfmt2(:,:)=0.d0
        do ist=1,nstfv
          i=i+1
          z1=evecsv(i,j)
          if (spinsprl.and.ssdph) z1=z1*zq(ispn)
          if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
            if (.not.done(ist,jspn)) then
              call wavefmt(lradstp,lmaxvr,ias,ngk(jspn,ik), &
               apwalm(:,:,:,:,jspn),evecfv(:,ist,jspn),lmmaxvr, &
               wfmt1(:,:,ist,jspn))
              done(ist,jspn)=.true.
            end if
! add to spinor wavefunction
            call zfmtadd(nrc,nrci,z1,wfmt1(:,:,ist,jspn),wfmt2)
          end if
        end do
! convert to spherical coordinates
        call zbsht(nrc,nrci,wfmt2,wfmt3(:,:,ispn))
      end do
    else
! spin-unpolarised wavefunction
      call wavefmt(lradstp,lmaxvr,ias,ngk(1,ik),apwalm,evecfv(:,j,1),lmmaxvr, &
       wfmt2)
! convert to spherical coordinates
      call zbsht(nrc,nrci,wfmt2,wfmt3)
    end if
! add to density and magnetisation
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
            z1=wfmt3(itp,irc,1)
            z2=wfmt3(itp,irc,2)
            t1=dble(z1)**2+aimag(z1)**2
            t2=dble(z2)**2+aimag(z2)**2
            z1=conjg(z1)*z2
            rhomt(itp,ir,ias)=rhomt(itp,ir,ias)+t0*(t1+t2)
            magmt(itp,ir,ias,1)=magmt(itp,ir,ias,1)+t4*dble(z1)
            magmt(itp,ir,ias,2)=magmt(itp,ir,ias,2)+t4*aimag(z1)
            magmt(itp,ir,ias,3)=magmt(itp,ir,ias,3)+t0*(t1-t2)
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
            t1=dble(wfmt3(itp,irc,1))**2+aimag(wfmt3(itp,irc,1))**2
            t2=dble(wfmt3(itp,irc,2))**2+aimag(wfmt3(itp,irc,2))**2
            rhomt(itp,ir,ias)=rhomt(itp,ir,ias)+t0*(t1+t2)
            magmt(itp,ir,ias,1)=magmt(itp,ir,ias,1)+t0*(t1-t2)
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
        rhomt(1:lmmax,ir,ias)=rhomt(1:lmmax,ir,ias) &
         +t0*(dble(wfmt3(1:lmmax,irc,1))**2+aimag(wfmt3(1:lmmax,irc,1))**2)
        if (irc.eq.nrci) lmmax=lmmaxvr
      end do
    end if
!$OMP END CRITICAL
  end do
! end loop over atoms
end do
if (tevecsv) deallocate(wfmt1)
deallocate(apwalm,wfmt2,wfmt3)
!------------------------------------------------!
!     interstitial density and magnetisation     !
!------------------------------------------------!
allocate(wfir(ngtot,nspinor))
do j=1,nstsv
  if (abs(occsv(j,ik)).lt.epsocc) cycle
  t0=wkpt(ik)*occsv(j,ik)
  t3=t0/omega
  t4=2.d0*t3
  wfir(:,:)=0.d0
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    i=0
    do ispn=1,nspinor
      jspn=jspnfv(ispn)
      do ist=1,nstfv
        i=i+1
        z1=evecsv(i,j)
        if (abs(dble(z1))+abs(aimag(z1)).gt.epsocc) then
          do igk=1,ngk(jspn,ik)
            ifg=igfft(igkig(igk,jspn,ik))
            wfir(ifg,ispn)=wfir(ifg,ispn)+z1*evecfv(igk,ist,jspn)
          end do
        end if
      end do
    end do
  else
! spin-unpolarised wavefunction
    do igk=1,ngk(1,ik)
      ifg=igfft(igkig(igk,1,ik))
      wfir(ifg,1)=evecfv(igk,j,1)
    end do
  end if
! Fourier transform wavefunction to real-space
  do ispn=1,nspinor
    call zfftifc(3,ngridg,1,wfir(:,ispn))
  end do
! add to density and magnetisation
!$OMP CRITICAL
  if (spinpol) then
! spin-polarised
    if (ncmag) then
! non-collinear
      do ir=1,ngtot
        z1=wfir(ir,1)
        z2=wfir(ir,2)
        t1=dble(z1)**2+aimag(z1)**2
        t2=dble(z2)**2+aimag(z2)**2
        z1=conjg(z1)*z2
        rhoir(ir)=rhoir(ir)+t3*(t1+t2)
        magir(ir,1)=magir(ir,1)+t4*dble(z1)
        magir(ir,2)=magir(ir,2)+t4*aimag(z1)
        magir(ir,3)=magir(ir,3)+t3*(t1-t2)
      end do
    else
! collinear
      do ir=1,ngtot
        t1=dble(wfir(ir,1))**2+aimag(wfir(ir,1))**2
        t2=dble(wfir(ir,2))**2+aimag(wfir(ir,2))**2
        rhoir(ir)=rhoir(ir)+t3*(t1+t2)
        magir(ir,1)=magir(ir,1)+t3*(t1-t2)
      end do
    end if
  else
! spin-unpolarised
    rhoir(:)=rhoir(:)+t3*(dble(wfir(:,1))**2+aimag(wfir(:,1))**2)
  end if
!$OMP END CRITICAL
end do
deallocate(wfir)
call timesec(ts1)
!$OMP ATOMIC
timerho=timerho+ts1-ts0
return
end subroutine
!EOC
