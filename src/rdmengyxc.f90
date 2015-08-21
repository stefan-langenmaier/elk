
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rdmengyxc
! !INTERFACE:
subroutine rdmengyxc
! !USES:
use modrdm
use modmain
! !DESCRIPTION:
!   Calculates RDMFT exchange-correlation energy.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik3,ik1,ik2
integer ist2,ist1,iv(3)
real(8) t1,t2,t3
! allocatable arays
real(8), allocatable :: vnlijji(:,:,:)
! external functions
real(8) r3taxi,rfmtinp
external r3taxi,rfmtinp
! calculate the prefactor
if (rdmxctype.eq.0) then
  engyx=0.d0
  return
! Hartree-Fock functional
else if (rdmxctype.eq.1) then
  t1=0.5d0/occmax
! Power functional
else if (rdmxctype.eq.2) then
  if (spinpol) then
    t1=0.5d0
  else
    t1=(0.25d0)**rdmalpha
  end if
else
  write(*,*)
  write(*,'("Error(rdmengyxc): rdmxctype not defined : ",I8)') rdmxctype
  write(*,*)
  stop
end if
! exchange-correlation energy
engyx=0.d0
allocate(vnlijji(nstsv,nstsv,nkpt))
! start loop over non-reduced k-points
do ik1=1,nkptnr
  call rdmgetvnl_ijji(ik1,vnlijji)
! find the equivalent reduced k-point
  iv(:)=ivknr(:,ik1)
  ik2=ikmap(iv(1),iv(2),iv(3))
  do ist1=1,nstsv
! start loop over reduced k-points
   do ik3=1,nkpt
     do ist2=1,nstsv
! Hartree-Fock functional
        if (rdmxctype.eq.1) then
          t2=t1*wkpt(ik3)*occsv(ist2,ik3)*occsv(ist1,ik2)
! Power functional
        else if (rdmxctype.eq.2) then
          t3=occsv(ist2,ik3)*occsv(ist1,ik2)
          if ((ist2.eq.ist1).and. &
           (r3taxi(vkl(:,ik3),vklnr(:,ik1)).lt.epslat)) then
            t2=(0.5d0/occmax)*wkpt(ik3)*t3
          else
            if (t3.gt.0.d0) then
              t2=t1*wkpt(ik3)*t3**rdmalpha
            else
              t2=0.d0
            end if
          end if
        end if
        engyx=engyx-t2*vnlijji(ist2,ist1,ik3)
      end do
    end do
  end do
end do
deallocate(vnlijji)
return
end subroutine
!EOC
