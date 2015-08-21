
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symrfmt
! !INTERFACE:
subroutine symrfmt(lrstp,is,rot,rfmt,srfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp : radial step length (in,integer)
!   is    : species number (in,integer)
!   rot   : rotation matrix (in,real(3,3))
!   rfmt  : input muffin-tin function (in,real(lmmaxvr,nrmtmax))
!   srfmt : output muffin-tin function (out,real(lmmaxvr,nrmtmax))
! !DESCRIPTION:
!   Applies a symmetry operation (in the form of a rotation matrix) to a real
!   muffin-tin function. See the routine {\tt rotrflm}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!   Changed from rotzflm to rotrflm, December 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp
integer, intent(in) :: is
real(8), intent(in) :: rot(3,3)
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax)
real(8), intent(out) :: srfmt(lmmaxvr,nrmtmax)
! local variables
integer ir,irc,nri,nro,iro,ld
nri=0
irc=0
do ir=1,nrmt(is),lrstp
  irc=irc+1
  if (ir.le.nrmtinr(is)) then
    srfmt(lmmaxinr+1:,ir)=0.d0
    nri=irc
  end if
end do
! first point in the outer point of the muffin-tin
iro=lrstp*nri+1
! number of points in the outer part
nro=irc-nri
! rotate the function
ld=lmmaxvr*lrstp
call rotrflm(rot,lmaxinr,nri,ld,rfmt,srfmt)
call rotrflm(rot,lmaxvr,nro,ld,rfmt(:,iro),srfmt(:,iro))
return
end subroutine
!EOC

