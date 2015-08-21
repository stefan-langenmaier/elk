
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genwfsv
! !INTERFACE:
subroutine genwfsv(tsh,tgp,tocc,ngp,igpig,evalsvp,apwalm,evecfv,evecsv,wfmt, &
 ld,wfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh     : .true. if wfmt should be in spherical harmonic basis (in,logical)
!   tgp     : .true. if wfir should be in G+p-space, otherwise in real-space
!             (in,logical)
!   tocc    : .true. if only occupied wavefunctions are required (in,logical)
!   ngp     : number of G+p-vectors (in,integer(nspnfv))
!   igpig   : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
!   evalsvp : second-variational eigenvalue for every state (in,real(nstsv))
!   apwalm  : APW matching coefficients
!             (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
!   evecfv  : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv  : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   wfmt    : muffin-tin part of the wavefunctions for every state in spherical
!             coordinates (out,complex(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
!   ld      : leading dimension, at least ngkmax if tgp is .true. and ngrtot
!             if tgp is .false. (in,integer)
!   wfir    : interstitial part of the wavefunctions for every state
!             (out,complex(ld,nspinor,nstsv))
! !DESCRIPTION:
!   Calculates the second-variational spinor wavefunctions in both the
!   muffin-tin and interstitial regions for every state of a particular
!   $k$-point. The wavefunctions in both regions are stored on a real-space
!   grid. A coarse radial mesh is assumed in the muffin-tins with with angular
!   momentum cut-off of {\tt lmaxvr}. If {\tt tocc} is {\tt .true.}, then only
!   the occupied states (those below the Fermi energy) are calculated.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!   Updated for spin-spirals, June 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
logical, intent(in) :: tgp
logical, intent(in) :: tocc
integer, intent(in) :: ngp(nspnfv)
integer, intent(in) :: igpig(ngkmax,nspnfv)
real(8), intent(in) :: evalsvp(nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(8), intent(out) :: wfir(ld,nspinor,nstsv)
! local variables
integer ispn,jspn,is,ia,ias
integer ist,i,j,n,igp,ifg
real(8) t1
complex(8) zq(2),zt1
! automatic arrays
logical done(nstfv)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:)
complex(8), allocatable :: wfmt2(:,:)
!--------------------------------!
!     muffin-tin wavefunction    !
!--------------------------------!
if (tevecsv) allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv))
if (.not.tsh) allocate(wfmt2(lmmaxvr,nrcmtmax))
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
    done(:)=.false.
    do j=1,nstsv
      if ((.not.tocc).or.((tocc).and.(evalsvp(j).lt.efermi))) then
        if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
          wfmt(:,:,ias,:,j)=0.d0
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
                if (.not.done(ist)) then
                  if (tsh) then
! wavefunction returned in spherical harmonics
                    call wavefmt(lradstp,lmaxvr,is,ia,ngp(jspn), &
                     apwalm(:,:,:,:,jspn),evecfv(:,ist,jspn),lmmaxvr, &
                     wfmt1(:,:,ist))
                  else
                    call wavefmt(lradstp,lmaxvr,is,ia,ngp(jspn), &
                     apwalm(:,:,:,:,jspn),evecfv(:,ist,jspn),lmmaxvr,wfmt2)
! convert from spherical harmonics to spherical coordinates
                    call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtvr, &
                     lmmaxvr,wfmt2,lmmaxvr,zzero,wfmt1(:,:,ist),lmmaxvr)
                  end if
                  done(ist)=.true.
                end if
! add to spinor wavefunction
                call zaxpy(n,zt1,wfmt1(:,:,ist),1,wfmt(:,:,ias,ispn,j),1)
              end if
! end loop over first-variational states
            end do
! end loop over spin
          end do
        else
! spin-unpolarised wavefunction
          if (tsh) then
! wavefunction returned in spherical harmonics
            call wavefmt(lradstp,lmaxvr,is,ia,ngp,apwalm,evecfv(:,j,1), &
             lmmaxvr,wfmt(:,:,ias,1,j))
          else
            call wavefmt(lradstp,lmaxvr,is,ia,ngp,apwalm,evecfv(:,j,1), &
             lmmaxvr,wfmt2)
! convert from spherical harmonics to spherical coordinates
            call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtvr,lmmaxvr, &
             wfmt2,lmmaxvr,zzero,wfmt(:,:,ias,1,j),lmmaxvr)
          end if
        end if
      end if
! end loop over second-variational states
    end do
! end loops over atoms and species
  end do
end do
if (tevecsv) deallocate(wfmt1)
if (.not.tsh) deallocate(wfmt2)
!-----------------------------------!
!     interstitial wavefunction     !
!-----------------------------------!
t1=1.d0/sqrt(omega)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,ispn,jspn,ist) &
!$OMP PRIVATE(zt1,igp,ifg)
!$OMP DO
do j=1,nstsv
  wfir(:,:,j)=0.d0
  if ((.not.tocc).or.((tocc).and.(evalsvp(j).lt.efermi))) then
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
            if (tgp) then
! wavefunction in G+p-space
              do igp=1,ngp(jspn)
                wfir(igp,ispn,j)=wfir(igp,ispn,j)+zt1*evecfv(igp,ist,jspn)
              end do
            else
! wavefunction in real-space
              zt1=t1*zt1
              do igp=1,ngp(jspn)
                ifg=igfft(igpig(igp,jspn))
                wfir(ifg,ispn,j)=wfir(ifg,ispn,j)+zt1*evecfv(igp,ist,jspn)
              end do
            end if
          end if
        end do
      end do
    else
! spin-unpolarised wavefunction
      if (tgp) then
        do igp=1,ngp(1)
          wfir(igp,1,j)=evecfv(igp,j,1)
        end do
      else
        do igp=1,ngp(1)
          ifg=igfft(igpig(igp,1))
          wfir(ifg,1,j)=t1*evecfv(igp,j,1)
        end do
      end if
    end if
! Fourier transform wavefunction to real-space if required
    if (.not.tgp) then
      do ispn=1,nspinor
        call zfftifc(3,ngrid,1,wfir(:,ispn,j))
      end do
    end if
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC

