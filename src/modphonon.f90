
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
! number of vectors for writing out frequencies and eigenvectors
integer nphwrt
! vectors in lattice coordinates for writing out frequencies and eigenvectors
real(8), allocatable :: vqlwrt(:,:)
! Coulomb pseudopotential
real(8) mustar
! number of temperatures for the Eliashberg equations and thermal properties
integer ntemp

!------------------------------------!
!     phonon supercell variables     !
!------------------------------------!
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
! original inverse reciprocal lattice vector matrix
real(8) binv0(3,3)
! original number of atoms
integer natoms0(maxspecies)
integer natmtot0
! original atomic positions in Cartesian coordinates
real(8) atposc0(3,maxatoms,maxspecies)
! original G-vector grid sizes
integer ngrid0(3),ngrtot0
! G-vector classification index which allows the Hamiltonian to be rearranged in
! banded diagonal form
integer, allocatable :: igph(:)

!------------------------------------------!
!     phonon linear-response variables     !
!------------------------------------------!
! derivative of density w.r.t. atomic displacement
complex(8), allocatable :: drhomt(:,:,:),drhoir(:)
! derivative of magnetisation w.r.t. atomic displacement
complex(8), allocatable :: dmagmt(:,:,:,:),dmagir(:,:)
! derivative of the effective potential
complex(8), allocatable :: dveffmt(:,:,:),dveffir(:)
! plane wave part of effective potential derivative
complex(8), allocatable :: dveffpw(:)
! derivative of the exchange-correlation magnetic field
complex(8), allocatable :: dbxcmt(:,:,:,:),dbxcir(:,:)
! plane wave part of magnetic field derivative
complex(8), allocatable :: dbxcpw(:,:)
! k+q-vectors in lattice and Cartesian coordinates
real(8), allocatable :: vkql(:,:),vkqc(:,:)
! number of G+k+q-vectors for each k-point
integer, allocatable :: ngkq(:,:)
! index from G+k+q-vectors to G-vectors
integer, allocatable :: igkqig(:,:,:)
! G+k+q-vectors in lattice and Cartesian coordinates
real(8), allocatable :: vgkql(:,:,:,:),vgkqc(:,:,:,:)
! G+k+q-vector lengths
real(8), allocatable :: gkqc(:,:,:)
! (theta, phi) coordinates of G+k+q-vectors
real(8), allocatable :: tpgkqc(:,:,:,:)
! structure factors for the G+k+q-vectors
complex(8), allocatable :: sfacgkq(:,:,:,:)

end module
