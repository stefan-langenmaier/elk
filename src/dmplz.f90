
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dmplz(l,k1,p,r,w2,dmpol)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   l     : angular momentum (in,integer)
!   k1    : k-index of tensor moment (in,integer)
!   p     : p-index of tensor moment  (in,integer)
!   r     : r-index of tensor moment  (in,integer)
!   w2    : modulus square of k1-p-r tensmom (in,real)
!   dmpol : polarization (out,real)
! !DESCRIPTION:
!   Calculate the polarization of each tensor component of the density matrix,
!   see {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (F.Cricchio and L.Nordstrom)
!EOP
!BOC
implicit none
! input variables
integer,intent(in) :: l
integer, intent(in) :: k1
integer,intent(in) :: p
integer,intent(in) :: r
real(8),intent(in) :: w2
real(8),intent(out) :: dmpol
! local variables
integer g
real(8) nk1l,fact
! external functions
real(8) wigner3j,factr,factnm
external wigner3j,factr,factnm
g=k1+p+r
if (g.eq.0) then
  fact=sqrt(w2)
  dmpol=fact*(2*(2*l+1)-fact)
else
  if (mod(g,2).eq.0) then
    fact=wigner3j(k1,p,r,0,0,0)
  else
    fact=sqrt(factnm(g-2*p,1)*factnm(g-2*r,1)/ &
     factnm(g+1,g-2*k1))*factnm(g,2)/(factnm(g-2*k1,2)* &
     factnm(g-2*p,2)*factnm(g-2*r,2))
  end if
  nk1l=factnm(2*l,1)/sqrt(factnm(2*l-k1,1)*factnm(2*l+k1+1,1))
  dmpol=fact**2*(2*k1+1)*(2*r+1)*(2*l+1)*nk1l**2*w2
end if
return
end subroutine
!EOC

