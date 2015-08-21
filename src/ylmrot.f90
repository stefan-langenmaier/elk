
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: ylmrot
! !INTERFACE:
subroutine ylmrot(p,alpha,beta,gamma,lmax,ld,d)
! !INPUT/OUTPUT PARAMETERS:
!   p     : if p=-1 then the rotation matrix is improper (in,integer)
!   alpha : first Euler angle (in,real)
!   beta  : second Euler angle (in,real)
!   gamma : third Euler angle (in,real)
!   lmax  : maximum angular momentum (in,integer)
!   ld    : leading dimension (in,integer)
!   d     : complex spherical harmonic rotation matrix (out,complex(ld,*))
! !DESCRIPTION:
!   Returns the rotation matrix in the basis of complex spherical harmonics
!   given the three Euler angles, $(\alpha,\beta,\gamma)$, and the parity, $p$,
!   of the rotation. The matrix is given by the formula
!   $$ D^l_{m_1m_2}(\alpha,\beta,\gamma)=d^l_{m_1m_2}(\beta)
!    e^{-i(m_1\alpha+m_2\gamma)}, $$
!   where $d$ is the rotation matrix about the $y$-axis. For improper rotations,
!   i.e. those which are a combination of a rotation and inversion, $D$ is
!   modified with $D^l_{m_1m_2}\rightarrow(-1)^l D^l_{m_1m_2}$. See the routines
!   {\tt euler} and {\tt ylmroty}.
!
! !REVISION HISTORY:
!   Created December 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: p
real(8), intent(in) :: alpha
real(8), intent(in) :: beta
real(8), intent(in) :: gamma
integer, intent(in) :: lmax
integer, intent(in) :: ld
complex(8), intent(out) :: d(ld,*)
! local variables
integer lmmax,l,m1,m2,lm1,lm2
real(8) t1
! allocatable arrays
real(8), allocatable :: dy(:,:)
if (lmax.lt.0) then
  write(*,*)
  write(*,'("Error(ylmrot): lmax < 0 : ",I8)') lmax
  write(*,*)
  stop
end if
lmmax=(lmax+1)**2
allocate(dy(lmmax,lmmax))
! generate the rotation matrix about the y-axis
call ylmroty(beta,lmax,lmmax,dy)
! rotation by alpha and gamma, as well as inversion if required
lm1=0
do l=0,lmax
  do m1=-l,l
    lm1=lm1+1
    lm2=l**2
    do m2=-l,l
      lm2=lm2+1
      t1=-dble(m1)*alpha-dble(m2)*gamma
      d(lm1,lm2)=dy(lm1,lm2)*cmplx(cos(t1),sin(t1),8)
      if ((p.eq.-1).and.(mod(l,2).ne.0)) d(lm1,lm2)=-d(lm1,lm2)
    end do
  end do
end do
deallocate(dy)
return
end subroutine
!EOC

