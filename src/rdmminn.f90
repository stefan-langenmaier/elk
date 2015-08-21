
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmminn
! !INTERFACE:
subroutine rdmminn
! !USES:
use modrdm
use modmain
! !DESCRIPTION:
!   Minimizes the total energy w.r.t. occupation numbers. The steepest-descent
!   algorithm is used.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! allocatable arrays
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
integer ik,it,idm
real(8) ep,de
! parameter to check energy convergence
real(8), parameter :: eps=1.d-8
ep=0.d0
! generate and write non-local matrix elements
if (wrtvnlijji) call rdmputvnl_ijji
open(61,file='RDMN_ENERGY.OUT',action='WRITE',form='FORMATTED')
if (spinpol) then
  open(62,file='RDMN_MOMENT.OUT',action='WRITE',form='FORMATTED')
end if
! begin iteration loop
do it=1,maxitn
  write(*,'("Info(rdmminn): iteration ",I4," of ",I4)') it,maxitn
! get the new occupation numbers
  call rdmvaryn
! zero the density
  rhomt(:,:,:)=0.d0
  rhoir(:)=0.d0
! zero the magnetisation
  if (spinpol) then
    magmt(:,:,:,:)=0.d0
    magir(:,:)=0.d0
  end if
! compute the charge density and magnetisation with the new occupancies
  allocate(evecfv(nmatmax,nstfv))
  allocate(evecsv(nstsv,nstsv))
  do ik=1,nkpt
! get the eigenvectors from file
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
    write(62,'(I6,3G18.10)') it,momtot(1:ndmag)
    call flushifc(62)
  end if
! normalise the density
  call rhonorm
! calculate the Coulomb potential
  call potcoul
! calculate Coulomb potential matrix elements (RDM states)
  call genvmat(vclmt,vclir,vclmat)
! calculate the energy
  call rdmenergy
! check for convergence
  de=ep-engytot
  if (it.gt.1) then
    if (abs(de).lt.eps) goto 10
  end if
  ep=engytot
! write energy and convergence factor to a file
  write(61,'(I6,2G18.10)') it,engytot,de
  call flushifc(61)
! end iteration loop
end do
10 continue
close(61)
if (spinpol) close(62)
return
end subroutine
!EOC
