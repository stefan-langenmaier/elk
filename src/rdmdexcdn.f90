
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmdexcdn
! !INTERFACE:
subroutine rdmdexcdn(dedn)
! !USES:
use modrdm
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   dedn : energy derivative (inout,real(nstsv,nkpt))
! !DESCRIPTION:
!   Calculates the derivative of the exchange-correlation energy w.r.t.
!   occupation numbers and adds the result to the total.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
real(8), intent(inout) :: dedn(nstsv,nkpt)
! local variables
integer ik1,ik2,ik3
integer ist2,ist1,iv(3)
! parameter for calculating the functional derivatives
real(8), parameter :: eps=1.d-12
real(8) t1,t2,t3,t4
! allocatable arays
real(8), allocatable :: vnlijji(:,:,:)
! external functions
real(8) r3taxi
external r3taxi
if (rdmxctype.eq.0) return
! calculate the pre-factor
if (rdmxctype.eq.1) then
  t1=1.d0/occmax
! power functional
else if (rdmxctype.eq.2) then
  if (spinpol) then
    t1=rdmalpha
  else
    t1=2.d0*rdmalpha*(0.25d0)**rdmalpha
  end if
else
  write(*,*)
  write(*,'("Error(rdmdexcdn): rdmxctype not defined : ",I8)') rdmxctype
  write(*,*)
  stop
end if
allocate(vnlijji(nstsv,nstsv,nkpt))
! start loop over non-reduced k-points
do ik1=1,nkptnr
! get the non-local matrix
  call rdmgetvnl_ijji(ik1,vnlijji)
! find the equivalent reduced k-point
  iv(:)=ivknr(:,ik1)
  ik2=ikmap(iv(1),iv(2),iv(3))
  do ist1=1,nstsv
!  start loop over reduced k-points
    do ik3=1,nkpt
      do ist2=1,nstsv
! Hartree-Fock functional
        if (rdmxctype.eq.1) then
          t2=t1*occsv(ist1,ik2)
! power functional
        else if (rdmxctype.eq.2) then
          if ((ist2.eq.ist1).and. &
           (r3taxi(vkl(1,ik3),vklnr(1,ik1)).lt.epslat)) then
            t2=(1.d0/occmax)*occsv(ist1,ik2)
          else
            t3=max(occsv(ist2,ik3),eps)
            t4=max(occsv(ist1,ik2),eps)
            t2=t1*(t4**rdmalpha)/(t3**(1.d0-rdmalpha))
          end if
        end if
        dedn(ist2,ik3)=dedn(ist2,ik3)+t2*vnlijji(ist2,ist1,ik3)
      end do
    end do
  end do
end do
deallocate(vnlijji)
return
end subroutine
!EOC

