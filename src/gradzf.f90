
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradzf(nr,ld1,r,vgpc,ld2,zfmt,zfir,gzfmt,gzfir)
use modmain
implicit none
! arguments
integer, intent(in) :: nr(nspecies)
integer, intent(in) :: ld1
real(8), intent(in) :: r(ld1,nspecies)
real(8), intent(in) :: vgpc(3,ngvec)
integer, intent(in) :: ld2
complex(8), intent(in) :: zfmt(lmmaxvr,ld2,natmtot),zfir(ngtot)
complex(8), intent(out) :: gzfmt(lmmaxvr,ld2,natmtot,3),gzfir(ngtot,3)
! local variables
integer is,ias,ig,ifg,i
complex(8) z1
! allocatable arrays
complex(8), allocatable :: gzfmt1(:,:,:),zfft(:)
! muffin-tin gradient
allocate(gzfmt1(lmmaxvr,ld2,3))
do ias=1,natmtot
  is=idxis(ias)
  call gradzfmt(lmaxvr,nr(is),r(:,is),lmmaxvr,ld2,zfmt(:,:,ias),gzfmt1)
  do i=1,3
    gzfmt(:,1:nr(is),ias,i)=gzfmt1(:,1:nr(is),i)
  end do
end do
deallocate(gzfmt1)
! interstitial gradient
allocate(zfft(ngtot))
zfft(:)=zfir(:)
call zfftifc(3,ngridg,-1,zfft)
do i=1,3
  gzfir(:,i)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
    z1=zfft(ifg)
    gzfir(ifg,i)=vgpc(i,ig)*cmplx(-aimag(z1),dble(z1),8)
  end do
  call zfftifc(3,ngridg,1,gzfir(:,i))
end do
deallocate(zfft)
return
end subroutine

