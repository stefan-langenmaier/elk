
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmatlu
use modmain
use modldapu
implicit none
! local variables
integer ik,ispn,ist
integer is,ia,ias
real(8) t1
! allocatable arrays
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
! zero the LDA+U density matrix
dmatlu(:,:,:,:,:)=0.d0
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecfv,evecsv,apwalm,dmat) &
!$OMP PRIVATE(ispn,is,ia,ias,ist,t1)
!$OMP DO
do ik=1,nkpt
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(dmat(lmmaxlu,lmmaxlu,nspinor,nspinor,nstsv))
! get the eigenvectors and occupancies from file
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
  call getoccsv(vkl(:,ik),occsv(:,ik))
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! begin loop over atoms and species
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      call gendmat(.false.,.false.,0,lmaxlu,is,ia,ngk(:,ik),apwalm,evecfv, &
       evecsv,lmmaxlu,dmat)
      do ist=1,nstsv
        t1=wkpt(ik)*occsv(ist,ik)
!$OMP CRITICAL
        dmatlu(:,:,:,:,ias)=dmatlu(:,:,:,:,ias)+t1*dmat(:,:,:,:,ist)
!$OMP END CRITICAL
      end do
    end do
  end do
  deallocate(evecfv,evecsv,apwalm,dmat)
end do
!$OMP END DO
!$OMP END PARALLEL
! symmetrise the density matrix
call symdmat(lmaxlu,lmmaxlu,dmatlu)
return
end subroutine

