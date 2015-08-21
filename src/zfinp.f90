
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfinp
! !INTERFACE:
complex(8) function zfinp(tsh,zfmt1,zfir1,zfmt2,zfir2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh   : .true. if the muffin-tin functions are in spherical harmonics
!           (in,logical)
!   zfmt1 : first complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngtot))
!   zfmt2 : second complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngtot))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions over the entire unit
!   cell. The muffin-tin functions should be stored on the coarse radial grid
!   and have angular momentum cut-off {\tt lmaxvr}. In the intersitial region,
!   the integrand is multiplied with the characteristic function, to remove the
!   contribution from the muffin-tin. See routines {\tt zfmtinp} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
complex(8), intent(in) :: zfmt1(lmmaxvr,nrcmtmax,natmtot),zfir1(ngtot)
complex(8), intent(in) :: zfmt2(lmmaxvr,nrcmtmax,natmtot),zfir2(ngtot)
! local variables
integer ias,is,ir
! external functions
complex(8) zfmtinp
external zfmtinp
zfinp=0.d0
! interstitial contribution
do ir=1,ngtot
  zfinp=zfinp+cfunir(ir)*conjg(zfir1(ir))*zfir2(ir)
end do
zfinp=zfinp*(omega/dble(ngtot))
! muffin-tin contribution
do ias=1,natmtot
  is=idxis(ias)
  zfinp=zfinp+zfmtinp(tsh,nrcmt(is),nrcmtinr(is),rcmt(:,is),zfmt1(:,:,ias), &
   zfmt2(:,:,ias))
end do
return
end function
!EOC

