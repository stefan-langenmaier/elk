
! Copyright (C) 2007-2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genkmatc(tvclcr)
use modmain
use modmpi
implicit none
! arguments
logical, intent(in) :: tvclcr
! local variables
integer ik,ist,ispn
integer is,ias,n,lp
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vmt(:,:,:),vir(:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: a(:,:)
! allocate global kinetic matrix elements array
if (allocated(kmatc)) deallocate(kmatc)
allocate(kmatc(nstsv,nstsv,nkpt))
allocate(vmt(lmmaxvr,nrcmtmax,natmtot),vir(ngtot))
! convert muffin-tin Kohn-Sham potential to spherical coordinates
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call rbsht(nrcmt(is),nrcmtinr(is),lradstp,vsmt(:,:,ias),1,vmt(:,:,ias))
end do
!$OMP END DO
!$OMP END PARALLEL
! mulitply Kohn-Sham potential by characteristic function
vir(:)=vsir(:)*cfunir(:)
! generate the spin-orbit coupling radial functions
call gensocfr
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv) &
!$OMP PRIVATE(wfmt,wfir,a,ispn,ist)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
  allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir(ngtot,nspinor,nstsv))
  allocate(a(nstsv,nstsv))
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
  call genwfsv(.false.,.false.,nstsv,idx,ngk(:,ik),igkig(:,:,ik),apwalm, &
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
    call genbmatk(bsmt,bsir,wfmt,wfir,a)
    kmatc(:,:,ik)=kmatc(:,:,ik)-a(:,:)
  end if
! add the Coulomb core matrix elements if required
  if (tvclcr) then
    call vclcore(wfmt,a)
    kmatc(:,:,ik)=kmatc(:,:,ik)+a(:,:)
  end if
! rotate kinetic matrix elements to Cartesian basis
  call zgemm('N','C',nstsv,nstsv,nstsv,zone,kmatc(:,:,ik),nstsv,evecsv,nstsv, &
   zzero,a,nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,a,nstsv,zzero, &
   kmatc(:,:,ik),nstsv)
  deallocate(apwalm,evecfv,evecsv)
  deallocate(wfmt,wfir,a)
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

