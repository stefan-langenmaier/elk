
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: tensmomengy
! !INTERFACE:
subroutine tensmomengy(is,l,k1,p,r,edir,exch)
! !USES:
use modldapu
! !INPUT/OUTPUT PARAMETERS:
!   is   : species type (in,integer)
!   l    : angular momentum (in,integer)
!   k1   : k-index of tensor moment (in,integer)
!   p    : p-index of tensor moment  (in,integer)
!   r    : r-index of tensor moment  (in,integer)
!   edir : direct energy (out,real)
!   exch : exchange energy (out,real)
! !DESCRIPTION:
!   Decomposes the LDA+U Hartree-Fock energy, direct and exchange term, in
!   tensor moments, see {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (F.Cricchio and L.Nordstrom)
!EOP
!BOC
implicit none
integer, intent(in) :: is
integer, intent(in) :: l
integer, intent(in) :: k1
integer, intent(in) :: p
integer, intent(in) :: r
real(8), intent(out) :: edir
real(8), intent(out) :: exch
! local variables
integer k,g
real(8) nk1l,fact
! automatic arrays
real(8) ie(0:2*lmaxlu),be(0:2*lmaxlu)
! external functions
real(8) factr,factnm,wigner3j,wigner6j
external factr,factnm,wigner3j,wigner6j
! calculate tensor moment decomposition of direct and exchange energy
! see Eq. 23 and 29 Paper and Eq.3-4 PRB 78, 100404 (2008)
g=k1+p+r
if (mod(g,2).eq.0) then
  fact=wigner3j(k1,p,r,0,0,0)
else
  fact=sqrt(factnm(g-2*p,1)*factnm(g-2*r,1)/ &
   factnm(g+1,g-2*k1))*factnm(g,2)/(factnm(g-2*k1,2)* &
   factnm(g-2*p,2)*factnm(g-2*r,2))
end if
nk1l=factnm(2*l,1)/sqrt(factnm(2*l-k1,1)*factnm(2*l+k1+1,1))
fact=fact**2*(2*r+1)
do k=0,2*l,2
  ie(k)=0.5d0*((2*l+1)*nk1l*wigner3j(l,k,l,0,0,0))**2
  be(k)=0.5d0*((2*k1+1)*(-1)**k1*wigner6j(l,l,k1,l,l,k)*fact)
end do
edir=0.d0
exch=0.d0
do k=0,2*l,2
  exch=exch-ie(k)*be(k)*flu(k,is)
end do
if ((p.eq.0).and.(mod(k1,2).eq.0)) then
  edir=ie(k1)*flu(k1,is)
end if
return
end subroutine
!EOC

