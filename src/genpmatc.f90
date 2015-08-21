
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpmatc
! !INTERFACE:
subroutine genpmatc
! !USES:
use modmain
use modmpi
use modtddft
! !DESCRIPTION:
!   Generates the momentum matrix, $P_C$, in the `Cartesian' basis. The momentum
!   matrix elements in the second-variational basis can then be obtained by
!   $$ P_{ij}={\bf c}_i^{\dag}P_C {\bf c}_j, $$
!   where ${\bf c}_i$ is the $i$th second-variational eigenvector. See also
!   {\tt genpmat}.
!
! !REVISION HISTORY:
!   Created May 2012 (K. Krieger)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,ispn
integer lp,n,i
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: a(:,:),b(:,:)
! allocate global Cartesian momentum matrix elements array
if (allocated(pmatc)) deallocate(pmatc)
allocate(pmatc(3,nstsv,nstsv,nkpt))
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(apwalm,evecfv,evecsv) &
!$OMP PRIVATE(wfmt,wfir,a,ispn,i)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
  allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir(ngkmax,nspinor,nstsv),a(nstsv,nstsv),b(nstsv,nstsv))
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! generate the second-variational wavefunctions for all states
  call genwfsv(.true.,.true.,nstsv,idx,ngk(:,ik),igkig(:,:,ik),apwalm,evecfv, &
   evecsv,wfmt,ngkmax,wfir)
! generate moment matrix elements
  call genpmat(ngk(:,ik),igkig(:,:,ik),vgkc(:,:,:,ik),wfmt,wfir,pmatc(:,:,:,ik))
! rotate matrix elements to Cartesian basis
  do i=1,3
    a(:,:)=pmatc(i,:,:,ik)
    call zgemm('N','C',nstsv,nstsv,nstsv,zone,a,nstsv,evecsv,nstsv,zzero,b, &
     nstsv)
    call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,b,nstsv,zzero,a, &
     nstsv)
    pmatc(i,:,:,ik)=a(:,:)
  end do
  deallocate(apwalm,evecfv,evecsv)
  deallocate(wfmt,wfir,a,b)
end do
!$OMP END DO
!$OMP END PARALLEL
! broadcast matrix elements to every process
if (np_mpi.gt.1) then
  n=3*nstsv*nstsv
  do ik=1,nkpt
    lp=mod(ik-1,np_mpi)
    call mpi_bcast(pmatc(:,:,:,ik),n,mpi_double_complex,lp,mpi_comm_world, &
     ierror)
  end do
end if
return
end subroutine
!EOC

