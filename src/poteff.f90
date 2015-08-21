
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
integer is,ia,ias,ir
real(8) ts0,ts1
call timesec(ts0)
! compute the Coulomb potential
call potcoul
! compute the exchange-correlation potential
call potxc
! add Coulomb and exchange-correlation potentials together
! muffin-tin part
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmtinr(is)
      veffmt(1:lmmaxinr,ir,ias)=vclmt(1:lmmaxinr,ir,ias) &
       +vxcmt(1:lmmaxinr,ir,ias)
      veffmt(lmmaxinr:lmmaxvr,ir,ias)=0.d0
    end do
    do ir=nrmtinr(is)+1,nrmt(is)
      veffmt(:,ir,ias)=vclmt(:,ir,ias)+vxcmt(:,ir,ias)
    end do
  end do
end do
! interstitial part
veffir(:)=vclir(:)+vxcir(:)
call timesec(ts1)
timepot=timepot+ts1-ts0
return
end subroutine
!EOC
