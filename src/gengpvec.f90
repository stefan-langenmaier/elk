
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gengpvec
! !INTERFACE:
subroutine gengpvec(vpl,vpc,ngp,igpig,vgpl,vgpc)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   vpl   : p-point vector in lattice coordinates (in,real(3))
!   vpc   : p-point vector in Cartesian coordinates (in,real(3))
!   ngp   : number of G+p-vectors returned (out,integer)
!   igpig : index from G+p-vectors to G-vectors (out,integer(ngkmax))
!   vgpl  : G+p-vectors in lattice coordinates (out,real(3,ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (out,real(3,ngkmax))
! !DESCRIPTION:
!   Generates a set of ${\bf G+p}$-vectors for the input $p$-point with length
!   less than {\tt gkmax}. These are used as the plane waves in the APW
!   functions.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Removed spherical coordinate generation, May 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vpc(3)
integer, intent(out) :: ngp
integer, intent(out) :: igpig(ngkmax)
real(8), intent(out) :: vgpl(3,ngkmax)
real(8), intent(out) :: vgpc(3,ngkmax)
! local variables
integer ig,igp
real(8) v(3),t0,t1
t0=gkmax**2
igp=0
do ig=1,ngvec
  v(:)=vgc(:,ig)+vpc(:)
  t1=v(1)**2+v(2)**2+v(3)**2
  if (t1.lt.t0) then
    igp=igp+1
    if (igp.gt.ngkmax) then
      write(*,*)
      write(*,'("Error(gengpvec): number of G+p-vectors exceeds ngkmax")')
      write(*,*)
      stop
    end if
! index to G-vector
    igpig(igp)=ig
! G+p-vector in lattice coordinates
    vgpl(:,igp)=dble(ivg(:,ig))+vpl(:)
! G+p-vector in Cartesian coordinates
    vgpc(:,igp)=v(:)
  end if
end do
ngp=igp
return
end subroutine
!EOC
