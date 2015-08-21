
! Copyright (C) 2011 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genfxc(vchi0,eps0,eps,fxc)
use modmain
implicit none
! arguments
complex(8), intent(in) :: vchi0(ngrpa,ngrpa,nwrpa)
complex(8), intent(in) :: eps0(ngrpa,ngrpa,nwrpa)
complex(8), intent(in) :: eps(ngrpa,ngrpa,nwrpa)
complex(8), intent(out) :: fxc(ngrpa,ngrpa,nwrpa)
! local variables
integer ig,jg,iw
real(8) t1
! allocatable arrays
complex(8), allocatable :: a(:,:)
! allocate local arrays
allocate(a(ngrpa,ngrpa))
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
    fxc(ig,ig,:)=-(fxclrc(1)+fxclrc(2)*dble(wrpa(:))**2)/fourpi
  end do
case(2)
! bootstrap
  fxc(:,:,:)=0.d0
  t1=-1.d0/(dble(eps0(1,1,1))-1.d0)
  do ig=1,ngrpa
    do jg=1,ngrpa
      fxc(ig,jg,:)=t1*eps(ig,jg,1)
    end do
  end do
case default
  write(*,*)
  write(*,'("Error(genfxc): fxctype not defined : ",I8)') fxctype
  write(*,*)
  stop
end select
! pre- and post-multiply by v^1/2 chi0 v^1/2
do iw=1,nwrpa
  call zgemm('N','N',ngrpa,ngrpa,ngrpa,zone,vchi0(:,:,iw),ngrpa,fxc(:,:,iw), &
   ngrpa,zzero,a,ngrpa)
  call zgemm('N','N',ngrpa,ngrpa,ngrpa,zone,a,ngrpa,vchi0(:,:,iw),ngrpa,zzero, &
   fxc(:,:,iw),ngrpa)
end do
deallocate(a)
return
end subroutine

