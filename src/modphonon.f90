
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modphonon
use modmain

!--------------------------!
!     phonon variables     !
!--------------------------!
! current phonon q-point, species, atom and polarisation index
integer iqph,isph,iaph,ipph
! number of primitive unit cells in phonon supercell
integer nphsc
! Cartesian offset vectors for each primitive cell in the supercell
real(8) vphsc(3,maxatoms)
! phonon displacement distance
real(8) deltaph
! original lattice vectors
real(8) avec0(3,3)
! original inverse of lattice vector matrix
real(8) ainv0(3,3)
! original number of atoms
integer natoms0(maxspecies)
integer natmtot0
! original atomic positions in Cartesian coordinates
real(8) atposc0(3,maxatoms,maxspecies)
! original G-vector grid sizes
integer ngrid0(3)
integer ngrtot0
! original effective potentials
real(8), allocatable :: veffmt0(:,:,:)
real(8), allocatable :: veffir0(:)
! derivative of the density w.r.t. atomic displacement
complex(8), allocatable :: drhomt(:,:,:)
complex(8), allocatable :: drhoir(:)
! derivative of the magnetisation w.r.t. atomic displacement
complex(8), allocatable :: dmagmt(:,:,:,:)
complex(8), allocatable :: dmagir(:,:)
! derivative of the effective potential
complex(8), allocatable :: dveffmt(:,:,:)
complex(8), allocatable :: dveffir(:)
! plane wave part of effective potential derivative
complex(8), allocatable :: dveffpw(:)
! derivative of the exchange-correlation magnetic field
complex(8), allocatable :: dbxcmt(:,:,:,:)
complex(8), allocatable :: dbxcir(:,:)
! plane wave part of magnetic field derivative
complex(8), allocatable :: dbxcpw(:,:)
! number of vectors for writing out frequencies and eigenvectors
integer nphwrt
! vectors in lattice coordinates for writing out frequencies and eigenvectors
real(8), allocatable :: vqlwrt(:,:)
! Coulomb pseudopotential
real(8) mustar

end module
