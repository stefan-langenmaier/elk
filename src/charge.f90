
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: charge
! !INTERFACE:
subroutine charge
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total charges by integrating the
!   density.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir
real(8) t1
! automatic arrays
real(8) fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
! external functions
real(8) ddot
external ddto
! find the muffin-tin charges
chgmttot=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      fr(ir)=rhomt(1,ir,ias)*spr(ir,is)**2
    end do
    call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
    chgmt(ias)=fourpi*y00*gr(nrmt(is))
    chgmttot=chgmttot+chgmt(ias)
  end do
end do
! find the interstitial charge
t1=ddot(ngrtot,rhoir,1,cfunir,1)
chgir=t1*omega/dble(ngrtot)
! total calculated charge
chgcalc=chgmttot+chgir
! write calculated total charge to test file
call writetest(400,'calculated total charge',tol=1.d-6,rv=chgcalc)
return
end subroutine
!EOC
