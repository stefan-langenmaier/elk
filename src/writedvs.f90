
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writedvs
use modmain
use modphonon
use modstore
implicit none
! local variables
integer is
open(50,file='DVS'//trim(filext),action='WRITE',form='UNFORMATTED')
write(50) version
write(50) nspecies
write(50) lmmaxvr
do is=1,nspecies
  write(50) natoms0(is)
  write(50) nrcmt(is)
end do
write(50) ngridg0
write(50) dvsmt,dvsir
close(50)
return
end subroutine

