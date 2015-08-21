
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: olpistl
! !INTERFACE:
subroutine olpistl(tapp,ngp,igpig,v,o)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   v     : set of input vectors to which O is applied if tapp is .true.,
!           otherwise not referenced (in,complex(*))
!   o     : O applied to v if tapp is .true., otherwise it is the overlap matrix
!           (inout,complex(*))
! !DESCRIPTION:
!   Computes the interstitial contribution to the overlap matrix for the APW
!   basis functions. The overlap is given by
!   $$ O^{\rm I}({\bf G+k,G'+k})=\tilde{\Theta}({\bf G-G'}), $$
!   where $\tilde{\Theta}$ is the characteristic function. See routine
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
complex(8), intent(in) :: v(*)
complex(8), intent(inout) :: o(*)
! local variables
integer ld,iv(3),ig
integer ist,i,j,k,ki,kj
! allocatable arrays
complex(8), allocatable :: oj(:)
if (tapp) then
! apply the overlap operator to v
  allocate(oj(ngp))
  do j=1,ngp
    do i=1,j
      iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
      ig=ivgig(iv(1),iv(2),iv(3))
      oj(i)=cfunig(ig)
    end do
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,k,ki,kj) SHARED(j)
!$OMP DO
    do ist=1,nstfv
      k=(ist-1)*nmatmax
      kj=k+j
      do i=1,j-1
        ki=k+i
        o(ki)=o(ki)+oj(i)*v(kj)
        o(kj)=o(kj)+conjg(oj(i))*v(ki)
      end do
      o(kj)=o(kj)+dble(oj(j))*v(kj)
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do
  deallocate(oj)
else
! calculate the matrix elements
  ld=ngp+nlotot
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,k,iv,ig)
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
      o(k)=o(k)+cfunig(ig)
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL
end if
return
end subroutine
!EOC

