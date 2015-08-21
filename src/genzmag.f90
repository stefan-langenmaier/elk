
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzmag(wfmt1,wfmt2,wfir1,wfir2,zmagmt,zmagir)
use modmain
implicit none
! arguments
complex(8), intent(in) ::  wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfir1(ngrtot,nspinor)
complex(8), intent(in) ::  wfir2(ngrtot,nspinor)
complex(8), intent(out) :: zmagmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
complex(8), intent(out) :: zmagir(ngrtot,ndmag)
! local variables
integer is,ia,ias,ir
complex(8) zt1,zt2
!-------------------------!
!     muffin-tin part     !
!-------------------------!
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    call genzmagmt(is,wfmt1(:,:,ias,1),wfmt1(:,:,ias,2),wfmt2(:,:,ias,1), &
     wfmt2(:,:,ias,2),natmtot,zmagmt(:,:,ias,1))
  end do
end do
!---------------------------!
!     interstitial part     !
!---------------------------!
! calculate the z-component of magnetisation: up-up - dn-dn
zmagir(:,ndmag)=conjg(wfir1(:,1))*wfir2(:,1)-conjg(wfir1(:,2))*wfir2(:,2)
! non-collinear case
if (ncmag) then
  do ir=1,ngrtot
! up-dn spin density
    zt1=conjg(wfir1(ir,1))*wfir2(ir,2)
! dn-up spin density
    zt2=conjg(wfir1(ir,2))*wfir2(ir,1)
! calculate the x-component: up-dn + dn-up
    zmagir(ir,1)=zt1+zt2
! calculate the y-component: i*(dn-up - up-dn)
    zt1=zt2-zt1
    zmagir(ir,2)=cmplx(-aimag(zt1),dble(zt1),8)
  end do
end if
return
end subroutine

