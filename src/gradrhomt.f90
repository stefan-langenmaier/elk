
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gradrhomt
use modmain
use modphonon
implicit none
! local variables
integer nr
! allocatable arrays
real(8), allocatable :: rfmt(:,:),grfmt(:,:,:)
! add gradient contribution from rigid shift of muffin-tin
allocate(rfmt(lmmaxvr,nrmtmax),grfmt(lmmaxvr,nrmtmax,3))
nr=nrmt(isph)
call gradrfmt(lmaxvr,nr,spr(:,isph),lmmaxvr,nrmtmax,rhomt(:,:,iasph),grfmt)
! convert to spherical coordinates
call rbsht(nr,nrmtinr(isph),grfmt(:,:,ipph),rfmt)
! subtract from density derivative
drhomt(:,1:nr,iasph)=drhomt(:,1:nr,iasph)-rfmt(:,1:nr)
deallocate(rfmt,grfmt)
return
end subroutine

