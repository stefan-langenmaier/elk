
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ggamt_2a
! !INTERFACE:
subroutine ggamt_2a(ias,g2rho,gvrho,grho2)
! !USES:
use modmain
! !DESCRIPTION:
!   Spin-unpolarised version of {\tt ggamt\_sp\_2a}.
!
! !REVISION HISTORY:
!   Created November 2009 (JKD and TMcQ)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ias
real(8), intent(out) :: g2rho(lmmaxvr,nrmtmax)
real(8), intent(out) :: gvrho(lmmaxvr,nrmtmax,3)
real(8), intent(out) :: grho2(lmmaxvr,nrmtmax)
! local variables
integer is,nr,i
! allocatable arrays
real(8), allocatable :: rfmt(:,:),grfmt(:,:,:)
allocate(rfmt(lmmaxvr,nrmtmax),grfmt(lmmaxvr,nrmtmax,3))
is=idxis(ias)
nr=nrmt(is)
! compute grad^2 rho in spherical coordinates
call grad2rfmt(lmaxvr,nr,spr(:,is),lmmaxvr,rhomt(:,:,ias),rfmt)
call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rfmt,lmmaxvr,0.d0, &
 g2rho,lmmaxvr)
! compute grad rho in spherical coordinates
call gradrfmt(lmaxvr,nr,spr(:,is),lmmaxvr,nrmtmax,rhomt(:,:,ias),grfmt)
do i=1,3
  call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,grfmt(:,:,i), &
   lmmaxvr,0.d0,gvrho(:,:,i),lmmaxvr)
end do
! (grad rho)^2
grho2(:,1:nr)=gvrho(:,1:nr,1)**2+gvrho(:,1:nr,2)**2+gvrho(:,1:nr,3)**2
deallocate(rfmt,grfmt)
return
end subroutine
!EOC

