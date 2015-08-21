
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmdkdc
! !INTERFACE:
subroutine rdmdkdc
! !USES:
use modrdm
use modmain
! !DESCRIPTION:
!   Calculates the derivative of kinetic energy w.r.t. the second-variational
!   coefficients {\tt evecsv}.
!
! !REVISION HISTORY:
!   Created October 2008 (Sharma)
!EOP
!BOC
implicit none
! allocatable arrays
complex(8), allocatable :: evecsv(:,:)
integer ik
allocate(evecsv(nstsv,nstsv))
do ik=1,nkpt
  call getevecsv(vkl(:,ik),evecsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,kinmatc(:,:,ik),nstsv,evecsv, &
   nstsv,zzero,dkdc(:,:,ik),nstsv)
end do
deallocate(evecsv)
return
end subroutine
!EOC
