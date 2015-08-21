
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rfmtctof
! !INTERFACE:
subroutine rfmtctof(rfmt)
! !INPUT/OUTPUT PARAMETERS:
!   rfmt : real muffin-tin function (in,real(lmmaxvr,nrmtmax,natmtot))
! !DESCRIPTION:
!   Converts a real muffin-tin function from a coarse to a fine radial mesh by
!   using cubic spline interpolation. See {\tt rfinterp} and {\tt spline}.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
use modmain
implicit none
! arguments
real(8), intent(inout) :: rfmt(lmmaxvr,nrmtmax,natmtot)
! local variables
integer ld,is,ias,lm
ld=lmmaxvr*lradstp
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is,lm)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
! interpolate with a clamped spline
  do lm=1,lmmaxvr
    call rfinterp(nrcmt(is),rcmt(:,is),ld,rfmt(lm,1,ias),nrmt(is),spr(:,is), &
     lmmaxvr,rfmt(lm,1,ias))
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC

