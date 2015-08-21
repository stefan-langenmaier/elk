
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modldapu
use modmain

!-------------------------!
!     LDA+U variables     !
!-------------------------!
! type of LDA+U to use (0: none)
integer ldapu
! input type for LDA+U calculation (1:5)
integer inptypelu
! maximum angular momentum
integer, parameter :: lmaxlu=3
integer, parameter :: lmmaxlu=(lmaxlu+1)**2
! angular momentum for each species
integer llu(maxspecies)
! U and J values for each species
real(8) ujlu(2,maxspecies)
! interpolation constant alpha for each atom [PRB 67, 153106 (2003)]
real(8), allocatable :: alphalu(:)
! readalu is .true. if alphalu is to be read from file
logical readalu
! tmomlu is .true. if tensor moments are written out at every iteration
logical tmomlu
! LDA+U density matrix
complex(8), allocatable :: dmatlu(:,:,:,:,:)
! LDA+U potential matrix in (l,m) basis
complex(8), allocatable :: vmatlu(:,:,:,:,:)
! LDA+U energy for each atom
real(8), allocatable :: engyalu(:)
! energy from the LDA+U correction
real(8) engylu
! Slater parameters
real(8) flu(0:2*lmaxlu,maxspecies)
! Racah parameters
real(8) elu(0:lmaxlu,maxspecies)
! screening length of Yukawa potential to calculate Slater integrals
real(8) lambdalu(maxspecies)
! energies to calculate radial functions for Slater integrals
real(8), allocatable :: flue(:,:,:)
! radial functions to calculate Slater integrals
real(8), allocatable :: flufr(:,:,:,:,:)
! fixed value of U for which screening length has to be determined
real(8) ulufix(maxspecies)
! initial values of screening length if U is fixed
real(8) lambdalu0(maxspecies)

end module

