
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmdedn
! !INTERFACE:
subroutine rdmdedn(dedn)
! !USES:
use modrdm
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   dedn : free energy derivative (out,real(nstsv,nkpt))
! !DESCRIPTION:
!   Calculates the negative of the derivative of total free energy w.r.t.
!   occupation numbers.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
real(8), intent(out) :: dedn(nstsv,nkpt)
! allocatable
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: c(:,:)
integer ik,ist
allocate(evecsv(nstsv,nstsv),c(nstsv,nstsv))
dedn(:,:)=0.d0
do ik=1,nkpt
! get evecsv from a file
  call getevecsv(vkl(:,ik),evecsv)
! kinetic contribution
  call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,dkdc(:,:,ik),nstsv, &
   zzero,c,nstsv)
  do ist=1,nstsv
! include Coulomb contribution
    dedn(ist,ik)=dedn(ist,ik)-(dble(c(ist,ist))+dble(vclmat(ist,ist,ik)))
  end do
end do
deallocate(evecsv,c)
! add exchange correlation contribution
call rdmdexcdn(dedn)
! add entropic contribution if needed
if (rdmtemp.gt.0.d0) then
  call rdmdtsdn(dedn)
end if
return
end subroutine
!EOC
