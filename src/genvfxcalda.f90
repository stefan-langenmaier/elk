
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvfxcalda(gqc,vfxc)
use modmain
implicit none
! arguments
real(8), intent(in) :: gqc(ngrpa)
complex(8), intent(out) :: vfxc(ngrpa,ngrpa,nwrpa)
! local variables
integer ig,jg,kg,iv(3)
real(8) t0
complex(8) zt1
! allocatable arrays
real(8), allocatable :: fxcmt(:,:,:),fxcir(:)
complex(8), allocatable :: fxcg(:)
allocate(fxcmt(lmmaxvr,nrmtmax,natmtot),fxcir(ngrtot))
allocate(fxcg(ngvec))
! generate the kernel f_xc in real-space
call kerxc(fxcmt,fxcir)
! Fourier transform the kernel
call zftrf(ngvec,ivg,vgc,fxcmt,fxcir,fxcg)
t0=1.d0/fourpi
do ig=1,ngrpa
  do jg=1,ngrpa
    iv(:)=ivg(:,ig)-ivg(:,jg)
    if ((iv(1).ge.intgv(1,1)).and.(iv(1).le.intgv(1,2)).and. &
        (iv(2).ge.intgv(2,1)).and.(iv(2).le.intgv(2,2)).and. &
        (iv(3).ge.intgv(3,1)).and.(iv(3).le.intgv(3,2))) then
      kg=ivgig(iv(1),iv(2),iv(3))
      zt1=t0*fxcg(kg)*(gqc(ig)*gqc(jg))
      vfxc(ig,jg,:)=zt1
    end if
  end do
end do
deallocate(fxcmt,fxcir,fxcg)
return
end subroutine

