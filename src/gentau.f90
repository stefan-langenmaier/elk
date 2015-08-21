
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gentau
! !INTERFACE:
subroutine gentau(taumt,tauir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   taumt : muffin-tin kinetic energy density (in,real(lmmaxvr,nrmtmax,natmtot))
!   tauir : interstitial kinetic energy density (in,real(ngrtot))
! !DESCRIPTION:
!   Computes the kinetic energy density
!   $$ \tau({\bf r})=\sum_{i{\bf k}}n_{i{\bf k}}\nabla\psi_{i{\bf k}}^*({\bf r})
!    \cdot\nabla\psi_{i{\bf k}}({\bf r}), $$
!   which includes a contraction over spinor components. Note that this is
!   \underline{twice} the value of the usual definition of the kinetic energy
!   density. Actually implimented is the equivalent but more efficient formula
!   $$ \tau({\bf r})=\sum_{i{\bf k}}n_{i{\bf k}}\varepsilon_{i{\bf k}}
!    |\psi_{i{\bf k}}({\bf r})|^2-V_{\rm eff}({\bf r})\rho({\bf r})
!    +\tfrac{1}{4}\nabla^2\rho({\bf r}). $$
!   See {\it J. Phys.: Condens. Matter} {\bf 19}, 196208 (2007).
!
! !REVISION HISTORY:
!   Created October 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(out) :: taumt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(out) :: tauir(ngrtot)
! local variables
integer ik,ld,is,ias
integer ir,irc,ig,ifg
! allocatable arrays
real(8), allocatable :: rfmt1(:,:),rfmt2(:,:)
complex(8), allocatable :: zfft(:)
! zero the kinetic energy density
taumt(:,:,:)=0.d0
tauir(:)=0.d0
! if wavefunctions do not exist tau cannot be computed
if ((iscl.le.1).and.(.not.trdstate)) return
! add contributions from each k-point
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkpt
  call gentauk(ik,taumt,tauir)
end do
!$OMP END DO
!$OMP END PARALLEL
! add the core kinetic energy density
call gentaucr(taumt)
ld=lmmaxvr*lradstp
! convert muffin-tin kinetic energy density to spherical harmonics
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rfmt1,is,irc,ir)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt1(lmmaxvr,nrcmtmax))
  is=idxis(ias)
  irc=0
  do ir=1,nrmt(is),lradstp
    irc=irc+1
    rfmt1(:,irc)=taumt(:,ir,ias)
  end do
  call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rfshtvr,lmmaxvr,rfmt1, &
   lmmaxvr,0.d0,taumt(:,:,ias),ld)
  deallocate(rfmt1)
end do
!$OMP END DO
!$OMP END PARALLEL
! symmetrise the kinetic energy density
call symrf(lradstp,taumt,tauir)
! convert from a coarse to a fine radial mesh
call rfmtctof(taumt)
! add density Laplacian term
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rfmt1,is,ir)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt1(lmmaxvr,nrmtmax))
  is=idxis(ias)
  call grad2rfmt(lmaxvr,nrmt(is),spr(:,is),lmmaxvr,rhomt(:,:,ias),rfmt1)
  do ir=1,nrmt(is)
    taumt(:,ir,ias)=taumt(:,ir,ias)+0.25d0*rfmt1(:,ir)
  end do
  deallocate(rfmt1)
end do
!$OMP END DO
!$OMP END PARALLEL
! Fourier transform interstitial density to G-space
allocate(zfft(ngrtot))
zfft(:)=rhoir(:)
call zfftifc(3,ngrid,-1,zfft)
! apply laplacian
do ig=1,ngrtot
  ifg=igfft(ig)
  zfft(ifg)=-zfft(ifg)*gc(ig)**2
end do
! Fourier transform back to real-space
call zfftifc(3,ngrid,1,zfft)
tauir(:)=tauir(:)+0.25d0*dble(zfft(:))
deallocate(zfft)
! convert back to spherical coordinates
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rfmt1,is,ir)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt1(lmmaxvr,nrmtmax))
  is=idxis(ias)
  do ir=1,nrmt(is)
    rfmt1(:,ir)=taumt(:,ir,ias)
  end do
  call dgemm('N','N',lmmaxvr,nrmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr,rfmt1, &
   lmmaxvr,0.d0,taumt(:,:,ias),lmmaxvr)
  deallocate(rfmt1)
end do
!$OMP END DO
!$OMP END PARALLEL
! subtract the density-potential term
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rfmt1,rfmt2,is,ir)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt1(lmmaxvr,nrmtmax),rfmt2(lmmaxvr,nrmtmax))
  is=idxis(ias)
  call dgemm('N','N',lmmaxvr,nrmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
   rhomt(:,:,ias),lmmaxvr,0.d0,rfmt1,lmmaxvr)
  call dgemm('N','N',lmmaxvr,nrmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
   veffmt(:,:,ias),lmmaxvr,0.d0,rfmt2,lmmaxvr)
  do ir=1,nrmt(is)
    taumt(:,ir,ias)=taumt(:,ir,ias)-rfmt1(:,ir)*rfmt2(:,ir)
  end do
  deallocate(rfmt1,rfmt2)
end do
!$OMP END DO
!$OMP END PARALLEL
tauir(:)=tauir(:)-rhoir(:)*veffir(:)
! multiply tau by 2 to conform to definition and make sure it is positive
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmt(is)
    taumt(:,ir,ias)=2.d0*abs(taumt(:,ir,ias))
  end do
end do
tauir(:)=2.d0*abs(tauir(:))
return
end subroutine
!EOC

