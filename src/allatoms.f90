
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: allatoms
! !INTERFACE:
subroutine allatoms
! !USES:
use modmain
! !DESCRIPTION:
!   Solves the Kohn-Sham-Dirac equations for each atom type in the solid and
!   finds the self-consistent radial wavefunctions, eigenvalues, charge
!   densities and potentials. The atomic densities can then be used to
!   initialise the crystal densities, and the atomic self-consistent potentials
!   can be appended to the muffin-tin potentials to solve for the core states.
!   Note that, irrespective of the value of {\tt xctype}, exchange-correlation
!   functional type 3 is used. See also {\tt atoms}, {\tt rhoinit},
!   {\tt gencore} and {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Modified for GGA, June 2007 (JKD)
!EOP
!BOC
implicit none
integer xctype_(3),xcgrad_
integer is
! allocatable arrays
real(8), allocatable :: rwf(:,:,:)
! allocate global species charge density and potential arrays
if (allocated(sprho)) deallocate(sprho)
allocate(sprho(spnrmax,nspecies))
if (allocated(spvr)) deallocate(spvr)
allocate(spvr(spnrmax,nspecies))
! use LDA (xctype=3) to set up atomic densities
xctype_(1)=3; xctype_(2:3)=0
xcgrad_=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rwf)
!$OMP DO
do is=1,nspecies
  allocate(rwf(spnrmax,2,spnstmax))
  call atom(solsc,ptnucl,spzn(is),spnst(is),spn(:,is),spl(:,is),spk(:,is), &
   spocc(:,is),xctype_,xcgrad_,spnr(is),spr(:,is),speval(:,is),sprho(:,is), &
   spvr(:,is),rwf)
  deallocate(rwf)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC

