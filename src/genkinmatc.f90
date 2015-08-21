
! Copyright (C) 2007-2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genkinmatc
use modmain
use modmpi
implicit none
! local variables
integer ld,is,ias
integer ik,ispn,ist,n,lp
! allocatable arrays
real(8), allocatable :: rfmt(:,:,:)
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
complex(8), allocatable :: c(:,:)
complex(8), allocatable :: bmat(:,:)
! allocate global kinetic matrix elements array
if (allocated(kinmatc)) deallocate(kinmatc)
allocate(kinmatc(nstsv,nstsv,nkpt))
! convert muffin-tin effective potential to spherical coordinates
allocate(rfmt(lmmaxvr,nrcmtmax,natmtot))
ld=lmmaxvr*lradstp
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
   veffmt(:,:,ias),ld,0.d0,rfmt(:,:,ias),lmmaxvr)
end do
!$OMP END DO
!$OMP END PARALLEL
! generate muffin-tin effective magnetic fields and s.o. coupling functions
call genbeffmt
! loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv,apwalm) &
!$OMP PRIVATE(wfmt,wfir,c,bmat,ispn,ist)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(evalfv(nstfv,nspnfv))
  allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
! solve the first- and second-variational secular equations
  call seceqn(ik,evalfv,evecfv,evecsv)
! write the first variational eigenvalues/vectors to file (this ensures the
! phase in eigenvectors is the same for subsequent matrix element evaluations)
  call putevalfv(ik,evalfv)
  call putevecfv(ik,evecfv)
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir(ngrtot,nspinor,nstsv))
  allocate(c(nstsv,nstsv))
  if (spinpol) allocate(bmat(nstsv,nstsv))
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! calculate the wavefunctions for all states of the input k-point
  call genwfsv(.false.,.false.,.false.,ngk(1,ik),igkig(:,1,ik),evalsv(:,ik), &
   apwalm,evecfv,evecsv,wfmt,ngrtot,wfir)
! compute effective potential matrix elements
  call genvmatk(rfmt,veffir,wfmt,wfir,kinmatc(:,:,ik))
  kinmatc(:,:,ik)=-kinmatc(:,:,ik)
! add second-variational eigenvalues along the diagonal
  do ist=1,nstsv
    kinmatc(ist,ist,ik)=kinmatc(ist,ist,ik)+evalsv(ist,ik)
  end do
! compute the exchange-correlation magnetic field matrix elements
  if (spinpol) then
    call genbmatk(beffmt,bxcir,wfmt,wfir,bmat)
    kinmatc(:,:,ik)=kinmatc(:,:,ik)-bmat(:,:)
  end if
! rotate kinetic matrix elements to Cartesian basis
  call zgemm('N','C',nstsv,nstsv,nstsv,zone,kinmatc(:,:,ik),nstsv,evecsv, &
   nstsv,zzero,c,nstsv)
  call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,c,nstsv,zzero, &
   kinmatc(:,:,ik),nstsv)
  deallocate(evalfv,evecfv,evecsv)
  deallocate(apwalm,wfmt,wfir,c)
  if (spinpol) deallocate(bmat)
end do
!$OMP END DO
!$OMP END PARALLEL
! broadcast matrix elements to every process
if (np_mpi.gt.1) then
  n=nstsv*nstsv
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(kinmatc(:,:,ik),n,mpi_double_complex,lp,mpi_comm_world, &
     ierror)
  end do
end if
deallocate(rfmt)
return
end subroutine

