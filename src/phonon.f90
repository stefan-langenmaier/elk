
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phonon
use modmain
use modphonon
implicit none
! local variables
integer ik,ispn,igkq
real(8) vl(3),vc(3)
! allocatable arrays
complex(8), allocatable :: dwfpw(:,:,:)
! initialise universal variables
call init0
call init1
call init2
! allocate local arrays
allocate(dwfpw(ngkmax,nspinor,nstsv))
! read density and potential from file
call readstate
! read in the eigenvalues and occupancies
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! allocate global arrays
if (allocated(vkql)) deallocate(vkql)
allocate(vkql(3,nkptnr))
if (allocated(vkqc)) deallocate(vkqc)
allocate(vkqc(3,nkptnr))
if (allocated(ngkq)) deallocate(ngkq)
allocate(ngkq(nspnfv,nkptnr))
if (allocated(igkqig)) deallocate(igkqig)
allocate(igkqig(ngkmax,nspnfv,nkptnr))
if (allocated(vgkql)) deallocate(vgkql)
allocate(vgkql(3,ngkmax,nspnfv,nkptnr))
if (allocated(vgkqc)) deallocate(vgkqc)
allocate(vgkqc(3,ngkmax,nspnfv,nkptnr))
if (allocated(gkqc)) deallocate(gkqc)
allocate(gkqc(ngkmax,nspnfv,nkptnr))
if (allocated(tpgkqc)) deallocate(tpgkqc)
allocate(tpgkqc(2,ngkmax,nspnfv,nkptnr))
if (allocated(sfacgkq)) deallocate(sfacgkq)
allocate(sfacgkq(ngkmax,natmtot,nspnfv,nkptnr))
10 continue
call dyntask(80)
if (iqph.eq.0) return
write(*,'("Info(phonon): working on ",A)') 'DYN'//trim(filext)
! loop over non-reduced k-point set
do ik=1,nkptnr
! k+q-vectors in lattice and Cartesian coordinates
  vkql(:,ik)=vkl(:,ik)+vql(:,iqph)
  vkqc(:,ik)=vkc(:,ik)+vqc(:,iqph)
  do ispn=1,nspnfv
    vl(:)=vkql(:,ik)
    vc(:)=vkqc(:,ik)
! spin-spiral case
    if (spinsprl) then
      if (ispn.eq.1) then
        vl(:)=vl(:)+0.5d0*vqlss(:)
        vc(:)=vc(:)+0.5d0*vqcss(:)
      else
        vl(:)=vl(:)-0.5d0*vqlss(:)
        vc(:)=vc(:)-0.5d0*vqcss(:)
      end if
    end if
! generate the G+k+q-vectors
    call gengpvec(vl,vc,ngkq(ispn,ik),igkqig(:,ispn,ik),vgkql(:,:,ispn,ik), &
     vgkqc(:,:,ispn,ik))
! generate the spherical coordinates of the G+k+q-vectors
    do igkq=1,ngkq(ispn,ik)
      call sphcrd(vgkqc(:,igkq,ispn,ik),gkqc(igkq,ispn,ik), &
       tpgkqc(:,igkq,ispn,ik))
    end do
! generate structure factors for the G+k+q-vectors
    call gensfacgp(ngkq(ispn,ik),vgkqc(:,:,ispn,ik),ngkmax,sfacgkq(:,:,ispn,ik))
  end do
end do
! begin the self-consistent loop
do iscl=1,maxscl
! zero the density and magnetisation derivatives
  drhomt(:,:,:)=0.d0
  drhoir(:)=0.d0
  if (spinpol) then
    dmagmt(:,:,:,:)=0.d0
    dmagir(:,:)=0.d0
  end if
!******
! end the self-consistent loop
end do
goto 10
return
end subroutine

