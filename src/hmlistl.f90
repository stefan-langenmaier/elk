
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlistl
! !INTERFACE:
subroutine hmlistl(tapp,ngp,igpig,vgpc,v,h)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tapp  : .true. if the Hamiltonian is to be applied to the input vector,
!           .false. if the full matrix is to be calculated (in,logical)
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   v     : set of input vectors to which H is applied if tapp is .true.,
!           otherwise not referenced (in,complex(*))
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
logical, intent(in) :: tapp
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: v(*)
complex(8), intent(inout) :: h(*)
! local variables
integer ld,iv(3),ig
integer ist,i,j,k,ki,kj
real(8) t1
! allocatable arrays
complex(8), allocatable :: hj(:)
if (tapp) then
! apply the Hamiltonian operator to v
  allocate(hj(ngp))
  do j=1,ngp
    do i=1,j
      iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
      ig=ivgig(iv(1),iv(2),iv(3))
      t1=0.5d0*(vgpc(1,i)*vgpc(1,j)+vgpc(2,i)*vgpc(2,j)+vgpc(3,i)*vgpc(3,j))
      hj(i)=veffig(ig)+t1*cfunig(ig)
    end do
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,k,ki,kj) SHARED(j)
!$OMP DO
    do ist=1,nstfv
      k=(ist-1)*nmatmax
      kj=k+j
      do i=1,j-1
        ki=k+i
        h(ki)=h(ki)+hj(i)*v(kj)
        h(kj)=h(kj)+conjg(hj(i))*v(ki)
      end do
      h(kj)=h(kj)+dble(hj(j))*v(kj)
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do
  deallocate(hj)
else
! calculate the matrix elements
  ld=ngp+nlotot
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,k,iv,ig,t1)
!$OMP DO
  do j=1,ngp
    if (tpmat) then
      k=((j-1)*j)/2
    else
      k=(j-1)*ld
    end if
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
end if
return
end subroutine
!EOC
