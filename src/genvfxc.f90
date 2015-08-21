
! Copyright (C) 2011 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvfxc(gqc,vchi0,eps0,eps,vfxc)
use modmain
use modtddft
implicit none
! arguments
real(8), intent(in) :: gqc(ngrf)
complex(8), intent(in) :: vchi0(nwrf,ngrf,ngrf)
complex(8), intent(in) :: eps0(ngrf,ngrf,nwrf)
complex(8), intent(in) :: eps(ngrf,ngrf,nwrf)
complex(8), intent(out) :: vfxc(ngrf,ngrf,nwrf)
! local variables
integer ig,jg,iw
real(8) t1
! allocatable arrays
complex(8), allocatable :: a(:,:),b(:,:)
! compute v^(-1/2) f_xc v^(-1/2)
select case(fxctype(1))
case(0,1)
! RPA
  vfxc(:,:,:)=0.d0
  return
case(3)
! ALDA
  call genvfxcg(gqc,vfxc)
case(200)
! long-range contribution with dynamic correlations
  vfxc(:,:,:)=0.d0
  do ig=1,ngrf
    vfxc(ig,ig,:)=-(fxclrc(1)+fxclrc(2)*dble(wrf(:))**2)/fourpi
  end do
case(210)
! bootstrap
  vfxc(:,:,:)=0.d0
  t1=-1.d0/(dble(eps0(1,1,1))-1.d0)
  do ig=1,ngrf
    do jg=1,ngrf
      vfxc(ig,jg,:)=t1*eps(ig,jg,1)
    end do
  end do
case default
  write(*,*)
  write(*,'("Error(genvfxc): fxctype not defined : ",I8)') fxctype
  write(*,*)
  stop
end select
! right multiply by v^1/2 chi0 v^1/2
allocate(a(ngrf,ngrf),b(ngrf,ngrf))
do iw=1,nwrf
  a(:,:)=vfxc(:,:,iw)
  b(:,:)=vchi0(iw,:,:)
  call zgemm('N','N',ngrf,ngrf,ngrf,zone,a,ngrf,b,ngrf,zzero,vfxc(:,:,iw),ngrf)
end do
deallocate(a,b)
return
end subroutine

