
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmminc
! !INTERFACE:
subroutine rdmminc
! !USES:
use modrdm
use modmain
! !DESCRIPTION:
!   Minimizes the total energy w.r.t. the second-variational coefficients
!   {\tt evecsv}. The steepest-descent algorithm is used.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
integer it,ik,idm
! parameter to check energy convergence
real(8), parameter :: eps=1.d-10
! allocatable arrays
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
! allocate arrays
open(91,file='RDMC_ENERGY.OUT',action='WRITE',form='FORMATTED')
if (spinpol) then
  open(92,file='RDMC_MOMENT.OUT',action='WRITE',form='FORMATTED')
end if
write(*,*)
! begin iteration loop
do it=1,maxitc
  write(*,'("Info(rdmminc): iteration ",I4," of ",I4)') it,maxitc
! vary evecsv and orthogonalise it
  call rdmvaryc
! zero the density
  rhomt(:,:,:)=0.d0
  rhoir(:)=0.d0
! zero the magnetisation
  if (spinpol) then
    magmt(:,:,:,:)=0.d0
    magir(:,:)=0.d0
  end if
  allocate(evecfv(nmatmax,nstfv))
  allocate(evecsv(nstsv,nstsv))
  do ik=1,nkpt
! get the eigenvectors and values from file
    call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    call getevecsv(vkl(:,ik),evecsv)
! add to the density and magnetisation
    call rhomagk(ik,evecfv,evecsv)
  end do
  deallocate(evecfv,evecsv)
! convert muffin-tin density/magnetisation to spherical harmonics
  call rhomagsh
! symmetrise the density
  call symrf(lradstp,rhomt,rhoir)
! convert the muffin-tin density from coarse to a fine grid
  call rfmtctof(rhomt)
  if (spinpol) then
! symmetrise the magnetisation
    call symrvf(lradstp,magmt,magir)
! convert the magnetisation from a coarse to a fine radial mesh
    do idm=1,ndmag
      call rfmtctof(magmt(:,:,:,idm))
    end do
  end if
! add core density to the valence density
  call addrhocr
! calculate the charges
  call charge
! calculate the magnetic moment
  if (spinpol) then
    call moment
    write(92,'(I6,3G18.10)') it,momtot(1:ndmag)
    call flushifc(92)
  end if
! normalise the density
  call rhonorm
! calculate the Coulomb potential
  call potcoul
! calculate Coulomb matrix elements
  call genvmat(vclmt,vclir,vclmat)
! calculate derivative of kinetic energy w.r.t. evecsv
  call rdmdkdc
! calculate the energy
  call rdmenergy
! write energy to a file
  write(91,'(I6,2G18.10)') it,engytot
  call flushifc(91)
! end iteration loop
end do
close(91)
close(92)
return
end subroutine
!EOC
