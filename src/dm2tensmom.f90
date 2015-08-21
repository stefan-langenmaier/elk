
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! !ROUTINE: dm2tensmom
! !INTERFACE:
subroutine dm2tensmom(l,nsp,ldim,k1,p,r,dmcomb,tmom)
! !INPUT/OUTPUT PARAMETERS:
!   l      : angular momentum (in,integer)
!   nsp    : number of spinors  (in,integer)
!   ldim   : ldim=2*l (in,integer)
!   k1     : k-index of tensor moment (in,integer)
!   p      : p-index of tensor moment  (in,integer)
!   r      : r-index of tensor moment  (in,integer)
!   dmcomb : density matrix with combined indices (in,complex)
!   tmom   : tensor moment (out,complex)
! !DESCRIPTION:
!   Transform the density matrix to tensor moment representation, see
!   {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (F.Cricchio and L.Nordstrom)
!EOP
!BOC
implicit none
integer, intent(in) :: l
integer, intent(in) :: nsp
integer, intent(in) :: ldim
integer, intent(in) :: k1
integer, intent(in) :: p
integer, intent(in) :: r
complex(8),intent(in) :: dmcomb(ldim*nsp,ldim*nsp)
complex(8),intent(out) :: tmom(-ldim:ldim)
! local variables
integer g,t,x,y,ispn,jspn
integer m1,m2,lm1,lm2,i,j
real(8) nk1l,nps,t1
complex(8), parameter :: zi=(0.d0,1.d0)
complex(8) zt1,zt2,zt3
! automatic arrays
real(8) w3js(0:1,-1:1,2,2)
! external functions
real(8) wigner3j,factr,factnm
external wigner3j,factr,factnm
! NOTE: wigner3j routine only handles integer indices
! tabulate the spin-3j symbols w3js(p,y,ispn,jspn) p=0,1, y=-p,p ispn,jspn=1,2
w3js(:,:,:,:)=0.d0
! 1 -> spin-up,   i.e. +1/2
! 2 -> spin-down, i.e. -1/2
w3js(0,0,1,1)= 1.d0/sqrt(2.d0)
w3js(0,0,2,2)=-1.d0/sqrt(2.d0)
w3js(1,0,1,1)= 1.d0/sqrt(6.d0)
w3js(1,0,2,2)= 1.d0/sqrt(6.d0)
w3js(1,-1,2,1)=-1.d0/sqrt(3.d0)
w3js(1,1,1,2)=-1.d0/sqrt(3.d0)
! calculate 3-indices tensor moment, see Eq. 23, 26, 27 in PRB 80, 035121 (2009)
nk1l=factnm(2*l,1)/sqrt(factnm(2*l-k1,1)*factnm(2*l+k1+1,1))
nps=1.d0/sqrt(factnm(2+p,1))
g=k1+p+r
if (mod(g,2).eq.0) then
  zt1=cmplx(1.d0/wigner3j(k1,p,r,0,0,0),0.d0,8)
else
  zt1=zi**(-g)*cmplx(sqrt(factr(g+1,g-2*k1)/factnm(g-2*p,1)/ &
   factnm(g-2*r,1)),0.d0,8)*cmplx(factnm(g-2*k1,2)* &
   factnm(g-2*p,2)*factnm(g-2*r,2)/factnm(g,2),0.d0,8)
end if
zt1=zt1*cmplx((-1)**(k1+p+l),0.d0,8)/cmplx(nk1l*nps,0.d0,8)
do t=-r,r
  tmom(t)=cmplx(0.d0,0.d0,8)
  do x=-k1,k1
    do y=-p,p
      zt2=cmplx(wigner3j(k1,r,p,-x,t,-y)*(-1)**(x+y),0.d0,8)
      do ispn=1,nsp
        do jspn=1,nsp
          zt3=cmplx(w3js(p,y,jspn,ispn),0.d0,8)
          if (abs(zt2).gt.1.d-10) then
            do m1=-l,l
              lm1=m1+l+1
              do m2=-l,l
                lm2=m2+l+1
                i=lm1+(2*l+1)*(ispn-1)
                j=lm2+(2*l+1)*(jspn-1)
                t1=wigner3j(l,k1,l,-m2,x,m1)
                tmom(t)=tmom(t)+zt1*zt2*zt3 &
                 *cmplx(t1*(-1)**(1+jspn-m2),0.d0,8)*dmcomb(i,j)
              end do
            end do
          end if
        end do
      end do
    end do
  end do
end do
return
end subroutine
!EOC

