
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!---------------------------------------------------------!
!     variables for storing original input parameters     !
!---------------------------------------------------------!

module modstore
use modmain

! lattice vectors
real(8) avec0(3,3)
! inverse reciprocal lattice vector matrix
real(8) binv0(3,3)
! number of atoms
integer natoms0(maxspecies)
integer natmtot0
! atomic positions in lattice coordinates
real(8) atposl0(3,maxatoms,maxspecies)
! atomic positions in Cartesian coordinates
real(8) atposc0(3,maxatoms,maxspecies)
! G-vector grid sizes
integer ngridg0(3),ngtot0
! external magnetic field in each muffin-tin
real(8) bfcmt00(3,maxatoms,maxspecies)
! muffin-tin fixed spin moments
real(8) mommtfix0(3,maxatoms,maxspecies)
! basis shifting flag
logical tshift0
! force calculation flag
logical tforce0
! automatic k-point generation flag
logical autokpt0
! primitive cell reduction flag
logical primcell0
! k-point grid
integer ngridk0(3)

end module

