
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

real(8) function rfint(rfmt,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: rfir(ngtot)
! local variables
integer is,ias,ir
real(8) sum
! automatic arrays
real(8) fr(nrmtmax),gr(nrmtmax)
sum=0.d0
! interstitial contribution
do ir=1,ngtot
  sum=sum+rfir(ir)*cfunir(ir)
end do
sum=sum*omega/dble(ngtot)
! muffin-tin contribution
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmt(is)
    fr(ir)=rfmt(1,ir,ias)*spr(ir,is)**2
  end do
  call fderiv(-1,nrmt(is),spr(:,is),fr,gr)
  sum=sum+fourpi*y00*gr(nrmt(is))
end do
rfint=sum
return
end function

