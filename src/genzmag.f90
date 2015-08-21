
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzmag(wfmt1,wfmt2,wfir1,wfir2,zmagmt,zmagir)
use modmain
implicit none
! arguments
complex(8), intent(in) ::  wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfir1(ngtot,nspinor)
complex(8), intent(in) ::  wfir2(ngtot,nspinor)
complex(8), intent(out) :: zmagmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
complex(8), intent(out) :: zmagir(ngtot,ndmag)
! local variables
integer is,ias,ir
complex(8) z1,z2
!-------------------------!
!     muffin-tin part     !
!-------------------------!
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call genzmagmt(is,wfmt1(:,:,ias,1),wfmt1(:,:,ias,2),wfmt2(:,:,ias,1), &
   wfmt2(:,:,ias,2),natmtot,zmagmt(:,:,ias,1))
end do
!$OMP END DO
!$OMP END PARALLEL
!---------------------------!
!     interstitial part     !
!---------------------------!
! calculate the z-component of magnetisation: up-up - dn-dn
zmagir(:,ndmag)=conjg(wfir1(:,1))*wfir2(:,1)-conjg(wfir1(:,2))*wfir2(:,2)
! non-collinear case
if (ncmag) then
  do ir=1,ngtot
! up-dn spin density
    z1=conjg(wfir1(ir,1))*wfir2(ir,2)
! dn-up spin density
    z2=conjg(wfir1(ir,2))*wfir2(ir,1)
! calculate the x-component: up-dn + dn-up
    zmagir(ir,1)=z1+z2
! calculate the y-component: i*(dn-up - up-dn)
    z1=z2-z1
    zmagir(ir,2)=cmplx(-aimag(z1),dble(z1),8)
  end do
end if
return
end subroutine

