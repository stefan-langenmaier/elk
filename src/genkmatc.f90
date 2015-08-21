
! Copyright (C) 2007-2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genkmatc(tvnlcr)
use modmain
use modmpi
implicit none
! arguments
logical, intent(in) :: tvnlcr
! local variables
integer ld,is,ias
integer ik,ispn,ist,n,lp
! allocatable arrays
real(8), allocatable :: vmt(:,:,:),vir(:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: c(:,:)
! allocate global kinetic matrix elements array
if (allocated(kmatc)) deallocate(kmatc)
allocate(kmatc(nstsv,nstsv,nkpt))
allocate(vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot))
! convert muffin-tin Kohn-Sham potential to spherical coordinates
ld=lmmaxvr*lradstp
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
   vsmt(:,:,ias),ld,0.d0,vmt(:,:,ias),lmmaxvr)
end do
!$OMP END DO
!$OMP END PARALLEL
! mulitply Kohn-Sham potential by characteristic function
vir(:)=vsir(:)*cfunir(:)
! generate the spin-orbit coupling radial functions
call gensocfr
! loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecfv,evecsv,apwalm) &
!$OMP PRIVATE(wfmt,wfir,c,ispn,ist)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir(ngtot,nspinor,nstsv))
  allocate(c(nstsv,nstsv))
! get the eigenvalues/vectors from file for input k-point
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! calculate the wavefunctions for all states of the input k-point
  call genwfsv(.false.,.false.,.false.,ngk(:,ik),igkig(:,:,ik),occsv,apwalm, &
   evecfv,evecsv,wfmt,ngtot,wfir)
! compute Kohn-Sham potential matrix elements
  call genvmatk(vmt,vir,wfmt,wfir,kmatc(:,:,ik))
  kmatc(:,:,ik)=-kmatc(:,:,ik)
! add second-variational eigenvalues along the diagonal
  do ist=1,nstsv
    kmatc(ist,ist,ik)=kmatc(ist,ist,ik)+evalsv(ist,ik)
  end do
! compute the exchange-correlation magnetic field matrix elements
  if (spinpol) then
    call genbmatk(bsmt,bsir,wfmt,wfir,c)
    kmatc(:,:,ik)=kmatc(:,:,ik)-c(:,:)
  end if
! add the non-local Coulomb core matrix elements if required
  if (tvnlcr) then
    call vnlcore(wfmt,c)
    kmatc(:,:,ik)=kmatc(:,:,ik)+c(:,:)
  end if
! rotate kinetic matrix elements to Cartesian basis
  call zgemm('N','C',nstsv,nstsv,nstsv,zone,kmatc(:,:,ik),nstsv,evecsv,nstsv, &
   zzero,c,nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,c,nstsv,zzero, &
   kmatc(:,:,ik),nstsv)
  deallocate(evecfv,evecsv)
  deallocate(apwalm,wfmt,wfir,c)
end do
!$OMP END DO
!$OMP END PARALLEL
! broadcast matrix elements to every process
if (np_mpi.gt.1) then
  n=nstsv*nstsv
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(kmatc(:,:,ik),n,mpi_double_complex,lp,mpi_comm_kpt,ierror)
  end do
end if
deallocate(vmt,vir)
return
end subroutine

