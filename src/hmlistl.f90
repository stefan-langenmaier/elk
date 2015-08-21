
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlistl
! !INTERFACE:
subroutine hmlistl(ngp,igpig,vgpc,h)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   h     : H applied to v if tapp is .true., otherwise it is the Hamiltonian
!           matrix (inout,complex(*))
! !DESCRIPTION:
!   Computes the interstitial contribution to the Hamiltonian matrix for the APW
!   basis functions. The Hamiltonian is given by
!   $$ H^{\rm I}({\bf G+k,G'+k})=\frac{1}{2}({\bf G+k})\cdot({\bf G'+k})
!    \tilde{\Theta}({\bf G-G'})+V^{\sigma}({\bf G-G'}), $$
!   where $V^{\sigma}$ is the effective interstitial potential for spin $\sigma$
!   and $\tilde{\Theta}$ is the characteristic function. See routine
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(inout) :: h(*)
! local variables
integer ld,iv(3),ig,i,j,k
real(8) t1
ld=ngp+nlotot
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,i,iv,ig,t1)
!$OMP DO
do j=1,ngp
  k=(j-1)*ld
  do i=1,j
    k=k+1
    iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
    ig=ivgig(iv(1),iv(2),iv(3))
    t1=0.5d0*(vgpc(1,i)*vgpc(1,j)+vgpc(2,i)*vgpc(2,j)+vgpc(3,i)*vgpc(3,j))
    h(k)=h(k)+veffig(ig)+t1*cfunig(ig)
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC
