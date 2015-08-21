
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: poteff
! !INTERFACE:
subroutine poteff
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the effective potential by adding together the Coulomb and
!   exchange-correlation potentials. See routines {\tt potcoul} and {\tt potxc}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,ir
real(8) ts0,ts1
call timesec(ts0)
! compute the Coulomb potential
call potcoul
! compute the exchange-correlation potential
call potxc
! effective potential from sum of Coulomb and exchange-correlation potentials
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmt(is)
    veffmt(:,ir,ias)=vclmt(:,ir,ias)+vxcmt(:,ir,ias)
  end do
end do
veffir(:)=vclir(:)+vxcir(:)
call timesec(ts1)
timepot=timepot+ts1-ts0
return
end subroutine
!EOC
