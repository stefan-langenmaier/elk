
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradrf(rfmt,rfir,grfmt,grfir)
use modmain
implicit none
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: rfir(ngrtot)
real(8), intent(out) :: grfmt(lmmaxvr,nrmtmax,natmtot,3)
real(8), intent(out) :: grfir(ngrtot,3)
! local variables
integer is,ias,i,ig,ifg
! allocatable arrays
real(8), allocatable :: grfmt1(:,:,:)
complex(8), allocatable :: zfft1(:),zfft2(:)
! muffin-tin gradient
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(grfmt1,is,i)
!$OMP DO
do ias=1,natmtot
  allocate(grfmt1(lmmaxvr,nrmtmax,3))
  is=idxis(ias)
  call gradrfmt(lmaxvr,nrmt(is),spr(:,is),lmmaxvr,nrmtmax,rfmt(:,:,ias),grfmt1)
  do i=1,3
    grfmt(:,1:nrmt(is),ias,i)=grfmt1(:,1:nrmt(is),i)
  end do
  deallocate(grfmt1)
end do
!$OMP END DO
!$OMP END PARALLEL
! interstitial gradient
allocate(zfft1(ngrtot),zfft2(ngrtot))
zfft1(:)=rfir(:)
call zfftifc(3,ngrid,-1,zfft1)
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    zfft2(ifg)=vgc(i,ig)*cmplx(-aimag(zfft1(ifg)),dble(zfft1(ifg)),8)
  end do
  call zfftifc(3,ngrid,1,zfft2)
  grfir(:,i)=dble(zfft2(:))
end do
deallocate(zfft1,zfft2)
return
end subroutine

