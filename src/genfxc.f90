
! Copyright (C) 2011 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genfxc(vchi0,fxc)
use modmain
implicit none
! arguments
complex(8), intent(in) :: vchi0(nwrpa,ngrpa,ngrpa)
complex(8), intent(out) :: fxc(nwrpa,ngrpa,ngrpa)
! local variables
integer ig,iw
! allocatable arrays
complex(8), allocatable :: a(:,:),b(:,:),c(:,:)
! allocate local arrays
allocate(a(ngrpa,ngrpa),b(ngrpa,ngrpa),c(ngrpa,ngrpa))
! compute v^(-1/2) fxc v^(-1/2)
select case(fxctype)
case(0)
! RPA
  fxc(:,:,:)=0.d0
  return
case(1)
! long-range contribution with dynamic correlations
  fxc(:,:,:)=0.d0
  do ig=1,ngrpa
    fxc(:,ig,ig)=-(fxclrc(1)+fxclrc(2)*dble(wrpa(:))**2)/fourpi
  end do
case default
  write(*,*)
  write(*,'("Error(genfxc): fxctype not defined : ",I8)') fxctype
  write(*,*)
  stop
end select
! pre- and post-multiply by v^1/2 chi0 v^1/2
do iw=1,nwrpa
  a(:,:)=vchi0(iw,:,:)
  b(:,:)=fxc(iw,:,:)
  call zgemm('N','N',ngrpa,ngrpa,ngrpa,zone,a,ngrpa,b,ngrpa,zzero,c,ngrpa)
  call zgemm('N','N',ngrpa,ngrpa,ngrpa,zone,c,ngrpa,a,ngrpa,zzero,b,ngrpa)
  fxc(iw,:,:)=b(:,:)
end do
deallocate(a,b,c)
return
end subroutine

