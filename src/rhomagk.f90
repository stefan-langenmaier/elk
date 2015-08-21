
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
!   {\tt rhoir}. See routines {\tt wavefmt}, {\tt genshtmat} and {\tt seceqn}.
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
integer ispn,jspn,ist
integer is,ia,ias,i,j,n
integer ir,irc,itp,igk,ifg
real(8) wo,t1,t2,t3
real(8) ts0,ts1
complex(8) zq(2),zt1,zt2,zt3
! automatic arrays
logical done(nstfv,nspnfv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt1(:,:)
complex(8), allocatable :: wfmt2(:,:,:,:)
complex(8), allocatable :: wfmt3(:,:,:)
complex(8), allocatable :: wfir(:,:)
call timesec(ts0)
!----------------------------------------------!
!     muffin-tin density and magnetisation     !
!----------------------------------------------!
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(wfmt1(lmmaxvr,nrcmtmax))
if (tevecsv) allocate(wfmt2(lmmaxvr,nrcmtmax,nstfv,nspnfv))
allocate(wfmt3(lmmaxvr,nrcmtmax,nspinor))
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
end do
do is=1,nspecies
  n=lmmaxvr*nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! de-phasing factor for spin-spirals
    if (spinsprl.and.ssdph) then
      t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
      zq(1)=cmplx(cos(t1),sin(t1),8)
      zq(2)=conjg(zq(1))
    end if
    done(:,:)=.false.
    do j=1,nstsv
      wo=wkpt(ik)*occsv(j,ik)
      if (abs(wo).gt.epsocc) then
        if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
          wfmt3(:,:,:)=0.d0
          i=0
          do ispn=1,nspinor
            if (spinsprl) then
              jspn=ispn
            else
              jspn=1
            end if
            do ist=1,nstfv
              i=i+1
              zt1=evecsv(i,j)
              if (spinsprl.and.ssdph) zt1=zt1*zq(ispn)
              if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
                if (.not.done(ist,jspn)) then
                  call wavefmt(lradstp,lmaxvr,is,ia,ngk(jspn,ik), &
                   apwalm(:,:,:,:,jspn),evecfv(:,ist,jspn),lmmaxvr,wfmt1)
! convert from spherical harmonics to spherical coordinates
                  call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtvr, &
                   lmmaxvr,wfmt1,lmmaxvr,zzero,wfmt2(:,:,ist,jspn),lmmaxvr)
                  done(ist,jspn)=.true.
                end if
! add to spinor wavefunction
                call zaxpy(n,zt1,wfmt2(:,:,ist,jspn),1,wfmt3(:,:,ispn),1)
              end if
            end do
          end do
        else
! spin-unpolarised wavefunction
          call wavefmt(lradstp,lmaxvr,is,ia,ngk(1,ik),apwalm,evecfv(:,j,1), &
           lmmaxvr,wfmt1)
! convert from spherical harmonics to spherical coordinates
          call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtvr,lmmaxvr, &
           wfmt1,lmmaxvr,zzero,wfmt3,lmmaxvr)
        end if
! add to density and magnetisation
!$OMP CRITICAL
        if (spinpol) then
! spin-polarised
          if (ncmag) then
! non-collinear
            irc=0
            do ir=1,nrmt(is),lradstp
              irc=irc+1
              do itp=1,lmmaxvr
                zt1=wfmt3(itp,irc,1)
                zt2=wfmt3(itp,irc,2)
                zt3=zt1*conjg(zt2)
                t1=dble(zt1)**2+aimag(zt1)**2
                t2=dble(zt2)**2+aimag(zt2)**2
                rhomt(itp,ir,ias)=rhomt(itp,ir,ias)+wo*(t1+t2)
                magmt(itp,ir,ias,1)=magmt(itp,ir,ias,1)+2.d0*wo*dble(zt3)
                magmt(itp,ir,ias,2)=magmt(itp,ir,ias,2)-2.d0*wo*aimag(zt3)
                magmt(itp,ir,ias,3)=magmt(itp,ir,ias,3)+wo*(t1-t2)
              end do
            end do
          else
! collinear
            irc=0
            do ir=1,nrmt(is),lradstp
              irc=irc+1
              do itp=1,lmmaxvr
                t1=dble(wfmt3(itp,irc,1))**2+aimag(wfmt3(itp,irc,1))**2
                t2=dble(wfmt3(itp,irc,2))**2+aimag(wfmt3(itp,irc,2))**2
                rhomt(itp,ir,ias)=rhomt(itp,ir,ias)+wo*(t1+t2)
                magmt(itp,ir,ias,1)=magmt(itp,ir,ias,1)+wo*(t1-t2)
              end do
            end do
          end if
        else
! spin-unpolarised
          irc=0
          do ir=1,nrmt(is),lradstp
            irc=irc+1
            rhomt(:,ir,ias)=rhomt(:,ir,ias) &
             +wo*(dble(wfmt3(:,irc,1))**2+aimag(wfmt3(:,irc,1))**2)
          end do
        end if
!$OMP END CRITICAL
      end if
    end do
  end do
end do
deallocate(apwalm,wfmt1,wfmt3)
if (tevecsv) deallocate(wfmt2)
!------------------------------------------------!
!     interstitial density and magnetisation     !
!------------------------------------------------!
allocate(wfir(ngrtot,nspinor))
do j=1,nstsv
  wo=wkpt(ik)*occsv(j,ik)
  if (abs(wo).gt.epsocc) then
    t3=wo/omega
    wfir(:,:)=0.d0
    if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
      i=0
      do ispn=1,nspinor
        if (spinsprl) then
          jspn=ispn
        else
          jspn=1
        end if
        do ist=1,nstfv
          i=i+1
          zt1=evecsv(i,j)
          if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
            do igk=1,ngk(jspn,ik)
              ifg=igfft(igkig(igk,jspn,ik))
              wfir(ifg,ispn)=wfir(ifg,ispn)+zt1*evecfv(igk,ist,jspn)
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
      call zfftifc(3,ngrid,1,wfir(:,ispn))
    end do
! add to density and magnetisation
!$OMP CRITICAL
    if (spinpol) then
! spin-polarised
      if (ncmag) then
! non-collinear
        do ir=1,ngrtot
          zt1=wfir(ir,1)
          zt2=wfir(ir,2)
          zt3=zt1*conjg(zt2)
          t1=dble(zt1)**2+aimag(zt1)**2
          t2=dble(zt2)**2+aimag(zt2)**2
          rhoir(ir)=rhoir(ir)+t3*(t1+t2)
          magir(ir,1)=magir(ir,1)+2.d0*t3*dble(zt3)
          magir(ir,2)=magir(ir,2)-2.d0*t3*aimag(zt3)
          magir(ir,3)=magir(ir,3)+t3*(t1-t2)
        end do
      else
! collinear
        do ir=1,ngrtot
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
  end if
end do
deallocate(wfir)
call timesec(ts1)
!$OMP CRITICAL
timerho=timerho+ts1-ts0
!$OMP END CRITICAL
return
end subroutine
!EOC
