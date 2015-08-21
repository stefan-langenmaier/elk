
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symmetry
use modmain
implicit none
! inverse of the lattice vector matrix
call r3minv(avec,ainv)
! find Bravais lattice symmetries
call findsymlat
! use only the identity if required
if (symtype.eq.0) nsymlat=1
! find the crystal symmetries and shift atomic positions if required
call findsymcrys
! find the site symmetries
call findsymsite
! check if fixed spin moments are invariant under the symmetry group
call checkfsm
! check if real symmetric eigen solver can be used
if (.not.tsyminv) tseqr=.false.
return
end subroutine

