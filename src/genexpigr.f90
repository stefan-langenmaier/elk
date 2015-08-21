
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genexpigr
use modmain
implicit none
! local variables
integer ig,ifg
if (allocated(expgmt)) deallocate(expgmt)
allocate(expgmt(lmmaxvr,nrcmtmax,natmtot,ngrpa))
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ig=1,ngrpa
  call genexpmt(vgc(:,ig),expgmt(:,:,:,ig))
end do
!$OMP END DO
!$OMP END PARALLEL
if (allocated(expgir)) deallocate(expgir)
allocate(expgir(ngrtot,ngrpa))
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ifg)
!$OMP DO
do ig=1,ngrpa
  ifg=igfft(ig)
  expgir(:,ig)=0.d0
  expgir(ifg,ig)=1.d0
  call zfftifc(3,ngrid,1,expgir(:,ig))
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

