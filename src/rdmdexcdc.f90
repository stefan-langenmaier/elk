
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmdexcdc
! !INTERFACE:
subroutine rdmdexcdc(dedc)
! !USES:
use modrdm
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   dedc : energy derivative (inout,complex(nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Calculates the derivative of the  exchange-correlation energy w.r.t.
!   {\tt evecsv} and adds the result to the total.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(inout) :: dedc(nstsv,nstsv,nkpt)
! local variables
integer ik1,ik2,jk,iv(3)
integer ist1,ist2,ist3,ist4
real(8) t1,t2
! allocatable arrays
complex(8), allocatable :: vnlijjk(:,:,:,:)
complex(8), allocatable :: evecsv(:,:)
! external functions
real(8) r3taxi
external r3taxi
if (rdmxctype.eq.0) return
! calculate the prefactor
if (rdmxctype.eq.1) then
! Hartree-Fock functional
  t1=1.d0/occmax
else if (rdmxctype.eq.2) then
! power functional
  if (spinpol) then
    t1=1.d0
  else
    t1=2.d0*(0.25d0)**rdmalpha
  end if
else
  write(*,*)
  write(*,'("Error(rdmdexcdc): rdmxctype not defined : ",I8)') rdmxctype
  write(*,*)
  stop
end if
allocate(vnlijjk(nstsv,nstsv,nstsv,nkpt))
allocate(evecsv(nstsv,nstsv))
! start loop over non-reduced k-points
do ik1=1,nkptnr
! get non-local matrix elements
  call rdmgetvnl_ijjk(ik1,vnlijjk)
! find the equivalent reduced k-point
  iv(:)=ivknr(:,ik1)
  jk=ikmap(iv(1),iv(2),iv(3))
  do ist4=1,nstsv
! start loop over reduced k-points
    do ik2=1,nkpt
! get the eigenvectors from file
      call getevecsv(vkl(:,ik2),evecsv)
      do ist3=1,nstsv
        do ist2=1,nstsv
          do ist1=1,nstsv
            if (rdmxctype.eq.1) then
! Hartree-Fock functional
              t2=t1*occsv(ist3,ik2)*occsv(ist4,jk)
            else if (rdmxctype.eq.2) then
! power functional
              if ((ist3.eq.ist4).and. &
               (r3taxi(vkl(:,ik2),vklnr(:,jk)).lt.epslat)) then
                t2=(1.d0/occmax)*occsv(ist4,jk)**2
              else
                t2=t1*(occsv(ist3,ik2)*occsv(ist4,jk))**rdmalpha
              end if
            end if
            dedc(ist2,ist3,ik2)=dedc(ist2,ist3,ik2)-t2*evecsv(ist2,ist1)* &
             vnlijjk(ist1,ist3,ist4,ik2)
          end do
        end do
      end do
! end loop over reduced k-points
    end do
  end do
! end loop over non-reduced k-points
end do
deallocate(vnlijjk,evecsv)
return
end subroutine
!EOC
