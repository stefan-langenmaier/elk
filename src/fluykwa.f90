
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: fluykwa
! !INTERFACE:
real(8) function fluykwa(is,l,k,lambda)
! !USES:
use modmain
use modldapu
! !INPUT/OUTPUT PARAMETERS:
!   is     : species type (in,integer)
!   l      : an angular momentum (in,integer)
!   k      : order of Slater parameter (in,integer)
!   lambda : screening length of Yukawa potential (in,real)
! !DESCRIPTION:
!   Calculates the Slater parameters using a screened Yukawa potential. See
!   {\it Phys. Rev. B} {\bf 52}, 1421 (1995) and {\it Phys. Rev. B} {\bf 80},
!   035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (Lars Nordstrom)
!   Modified and tested July 2008 (LN and FC)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: is
integer, intent(in) :: l
integer, intent(in) :: k
real(8), intent(in) :: lambda
! local variables
integer ias,io,nr,ir
integer nr1,nr2,ir1,ir2
real(8) r2,x,t1
! automatic arrays
real(8) clow(nrmtmax),chigh(nrmtmax)
real(8) blow(nrmtmax),bhigh(nrmtmax),fint(nrmtmax)
real(8) bint(nrmtmax),cint(nrmtmax),cf(4,nrmtmax)
! allocatable arrays
real(8), allocatable :: a(:,:)
real(8), allocatable :: b(:,:)
ias=idxas(1,is)
io=1
nr=nrmt(is)
! global factor Eq.5, Notes 30-05-08
! (-1)**k factor takes care of the additional (-1)**k introduced by zbessela(b)
t1=lambda*dble((2*k+1)*(-1)**(k+1))
! allocate Bessel and Hankel functions
allocate(a(0:2*l,nr),b(0:2*l,nr))
! zero all quantities
a(:,:)=0.d0
b(:,:)=0.d0
bint(:)=0.d0
blow(:)=0.d0
bhigh(:)=0.d0
clow(:)=0.d0
chigh(:)=0.d0
! calculate Slater parameters following Eq.6, Notes 29-04-09
do ir=1,nr
  r2=spr(ir,is)**2
  bint(ir)=flufr(ir,1,io,l,ias)*flufr(ir,1,io,l,ias)*r2
! argument of Bessel and Hankel functions divided by zi
  x=spr(ir,is)*lambda
! calculate Bessel and Hankel functions divided by zi
  call zbessela(2*l,x,a(0,ir))
  call zbesselb(2*l,x,b(0,ir))
end do
do ir=1,nr
  nr1=ir
  nr2=nr-ir+1
! 1st term: r1 < r
  do ir1=1,nr1
    ir2=ir1
    blow(ir1)=bint(ir2)*a(k,ir2)
  end do
! integrate 1st term
  call fderiv(-1,nr1,spr(1,is),blow,clow,cf)
! 2nd term : r2 > r
  do ir1=1,nr2
    ir2=ir1+ir-1
    bhigh(ir1)=bint(ir2)*b(k,ir2)
  end do
! integrate 2nd term
  call fderiv(-1,nr2,spr(ir,is),bhigh,chigh,cf)
! sum of the two terms
  cint(ir)=bint(ir)*(b(k,ir)*clow(nr1)+a(k,ir)*chigh(nr2))
end do
call fderiv(-1,nr,spr(1,is),cint,fint,cf)
fluykwa=t1*fint(nr)
deallocate(a,b)
return
end function
!EOC

