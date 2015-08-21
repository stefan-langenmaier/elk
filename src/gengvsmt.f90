
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengvsmt
use modmain
use modphonon
implicit none
! local variables
integer nr,ir
! allocatable arrays
complex(8), allocatable :: zfmt(:,:),gzfmt(:,:,:)
allocate(zfmt(lmmaxvr,nrmtmax),gzfmt(lmmaxvr,nrmtmax,3))
nr=nrmt(isph)
! convert potential to complex spherical harmonics
do ir=1,nr
  call rtozflm(lmaxvr,vsmt(:,ir,iasph),zfmt(:,ir))
end do
! calculate the gradient
call gradzfmt(lmaxvr,nr,spr(:,isph),lmmaxvr,nrmtmax,zfmt,gzfmt)
! copy current polarisation component to global array
gvsmt(:,1:nr)=gzfmt(:,1:nr,ipph)
deallocate(zfmt,gzfmt)
return
end subroutine
