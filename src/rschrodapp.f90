
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rschrodapp
! !INTERFACE:
subroutine rschrodapp(sol,l,nr,r,vr,p0,q0,q1,hp0)
! !INPUT/OUTPUT PARAMETERS:
!   sol : speed of light in atomic units (in,real)
!   l   : angular momentum quantum number (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   p0  : m th energy derivative of P (in,real(nr))
!   q0  : m th energy derivative of Q (in,real(nr))
!   q1  : radial derivative of q0 (in,real(nr))
!   hp0 : H applied to P (out,real(nr))
! !DESCRIPTION:
!   Applies the scalar relativistic radial Hamiltonian, $H$, to a radial
!   wavefunction, $P_l$. This is an approximation since we assume $P_l$ is a
!   scalar wavefunction, normalisable to unity. A Hamiltonian which satisfies
!   $H P_l=E P_l$ is given implicitly by
!   $$ H P_l=\left[\frac{l(l+1)}{2Mr^2}+V\right]P_l-\frac{1}{r}Q_l
!    -\frac{d}{dr}Q_l, $$
!   where $V$ is the external potential, $M=1-V/2c^2$ and $Q_l$ is obtained from
!   integrating the coupled scalar relativistic equations. See the routine
!   {\tt rschrodint} for further details.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
implicit none
real(8), intent(in) :: sol
integer, intent(in) :: l
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: vr(nr)
real(8), intent(in) :: p0(nr)
real(8), intent(in) :: q0(nr)
real(8), intent(in) :: q1(nr)
real(8), intent(out) :: hp0(nr)
! local variables
integer ir
real(8) t1,t2,t3
t1=1.d0/sol**2
t2=dble(l*(l+1))
do ir=1,nr
  t3=2.d0-t1*vr(ir)
  t3=t2/(t3*r(ir)**2)
  hp0(ir)=(t3+vr(ir))*p0(ir)-q0(ir)/r(ir)-q1(ir)
end do
return
end subroutine
!EOC
