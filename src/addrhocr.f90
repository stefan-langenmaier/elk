
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: addrhocr
! !INTERFACE:
subroutine addrhocr
! !USES:
use modmain
! !DESCRIPTION:
!   Adds the core density to the muffin-tin and interstitial densities. A
!   uniform background density is added in the interstitial region to take into
!   account leakage of core charge from the muffin-tin spheres.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias
integer nr,ir,ispn
real(8) sum1,sum2
real(8) v(3),t1,t2
! automatic arrays
real(8) fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
sum1=0.d0
sum2=0.d0
do is=1,nspecies
  nr=nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! loop over spin channels
    do ispn=1,nspncr
      do ir=1,nr
! add the core density to the muffin-tin density
        rhomt(1,ir,ias)=rhomt(1,ir,ias)+rhocr(ir,ias,ispn)/y00
        fr(ir)=rhocr(ir,ias,ispn)*spr(ir,is)**2
      end do
! compute the core charge inside the muffin-tins
      call fderiv(-1,nr,spr(:,is),fr,gr,cf)
      sum1=sum1+fourpi*gr(nr)
    end do
! add to the magnetisation in the case of a spin-polarised core
    if (spincore) then
      if (ncmag) then
! non-collinear
        do ir=1,nr
          t1=abs(rhocr(ir,ias,1)-rhocr(ir,ias,2))/y00
          v(:)=magmt(1,ir,ias,:)
          t2=sqrt(v(1)**2+v(2)**2+v(3)**2)
          if (t2.gt.1.d-10) magmt(1,ir,ias,:)=magmt(1,ir,ias,:)+(t1/t2)*v(:)
        end do
      else
! collinear
        do ir=1,nr
          t1=abs(rhocr(ir,ias,1)-rhocr(ir,ias,2))/y00
          if (magmt(1,ir,ias,1).lt.0.d0) t1=-t1
          magmt(1,ir,ias,1)=magmt(1,ir,ias,1)+t1
        end do
      end if
    end if
  end do
  sum2=sum2+dble(natoms(is))*(4.d0*pi/3.d0)*rmt(is)**3
end do
! add remaining core charge to interstitial density
chgcrlk=chgcr-sum1
t1=chgcrlk/(omega-sum2)
rhoir(:)=rhoir(:)+t1
return
end subroutine
!EOC

