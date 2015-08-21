
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zpotclmt
! !INTERFACE:
subroutine zpotclmt(ptnucl,lmax,nr,r,zn,ld,zrhomt,zvclmt)
! !INPUT/OUTPUT PARAMETERS:
!   ptnucl : .true. if the nucleus is a point particle (in,logical)
!   lmax   : maximum angular momentum (in,integer)
!   nr     : number of radial mesh points (in,integer)
!   r      : radial mesh (in,real(nr))
!   zn     : nuclear charge at the atomic center (in,real)
!   ld     : leading dimension (in,integer)
!   zrhomt : muffin-tin charge density (in,complex(ld,nr))
!   zvclmt : muffin-tin Coulomb potential (out,complex(ld,nr))
! !DESCRIPTION:
!   Solves the Poisson equation for the charge density contained in an isolated
!   muffin-tin using the Green's function approach. In other words, the
!   spherical harmonic expansion of the Coulomb potential, $V_{lm}$, is obtained
!   from the density expansion, $\rho_{lm}$, by
!   $$ V_{lm}(r)=\frac{4\pi}{2l+1}\left(\frac{1}{r^{l+1}}\int_0^r
!    \rho_{lm}(r'){r'}^{l+2}dr'+r^l\int_r^R\frac{\rho_{lm}(r')}{{r'}^{l-1}}dr'
!    \right)+\frac{1}{Y_{00}}\frac{z}{r}\delta_{l,0} $$
!   where the last term is the monopole arising from the point charge $z$, and
!   $R$ is the muffin-tin radius.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: ptnucl
integer, intent(in) :: lmax
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: zn
integer, intent(in) :: ld
complex(8), intent(in) :: zrhomt(ld,nr)
complex(8), intent(out) :: zvclmt(ld,nr)
! local variables
integer l,m,lm,ir
real(8), parameter :: fourpi=12.566370614359172954d0
! spherical harmonic for l=m=0
real(8), parameter :: y00=0.28209479177387814347d0
real(8) t1,t2,t3,t4
! automatic arrays
real(8) ri(nr),rl(4,nr),cf(4,nr)
real(8) fr1(nr),fr2(nr),fr3(nr),fr4(nr),fr5(nr)
! initialise r^l, r^(-l-1), r^(l+2) and r^(-l+1)
do ir=1,nr
  ri(ir)=1.d0/r(ir)
  rl(1,ir)=1.d0
  rl(2,ir)=ri(ir)
  t1=fourpi*r(ir)
  rl(3,ir)=t1*r(ir)
  rl(4,ir)=t1
end do
lm=0
do l=0,lmax
  do m=-l,l
    lm=lm+1
    do ir=1,nr
      t1=dble(zrhomt(lm,ir))
      t2=aimag(zrhomt(lm,ir))
      fr1(ir)=t1*rl(3,ir)
      fr2(ir)=t2*rl(3,ir)
      fr3(ir)=t1*rl(4,ir)
      fr4(ir)=t2*rl(4,ir)
    end do
    call fderiv(-1,nr,r,fr1,fr5,cf)
    call fderiv(-1,nr,r,fr2,fr1,cf)
    call fderiv(-1,nr,r,fr3,fr2,cf)
    call fderiv(-1,nr,r,fr4,fr3,cf)
    t1=fr2(nr)
    t2=fr3(nr)
    do ir=1,nr
      t3=rl(2,ir)*fr5(ir)+rl(1,ir)*(t1-fr2(ir))
      t4=rl(2,ir)*fr1(ir)+rl(1,ir)*(t2-fr3(ir))
      zvclmt(lm,ir)=cmplx(t3,t4,8)
    end do
  end do
  if (l.lt.lmax) then
    t1=fourpi/dble(2*(l+1)+1)
    do ir=1,nr
      rl(1,ir)=rl(1,ir)*r(ir)
      rl(2,ir)=rl(2,ir)*ri(ir)
      t2=t1*r(ir)**2
      rl(3,ir)=rl(1,ir)*t2
      rl(4,ir)=rl(2,ir)*t2
    end do
  end if
end do
! add the nuclear potential
if (zn.ne.0.d0) then
  call potnucl(ptnucl,nr,r,zn,fr1)
  t1=1.d0/y00
  zvclmt(1,:)=zvclmt(1,:)+t1*fr1(:)
end if
return
end subroutine
!EOC

