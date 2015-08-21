
! Copyright (C) 2008  F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeengyflu
! !INTERFACE:
subroutine writeengyflu
! !USES:
use modmain
use modldapu
! !DESCRIPTION:
!   Writes to file the linearisation energies for all radial functions used to
!   calculate the Slater integrals through a Yukawa potential.
!
! !REVISION HISTORY:
!   Created July 2008 (Francesco Cricchio)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,l,io
open(50,file='ENGYFLU'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  l=llu(is)
  io=1
  if (l.lt.0) goto 10
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    write(50,'(" Radial functions used for Slater parameters :")')
    write(50,'("  l = ",I2,", order = ",I2," : ",G18.10)') l,io,flue(io,l,ias)
  end do
10 continue
end do
close(50)
return
end subroutine
!EOC

