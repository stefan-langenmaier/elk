
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bse
use modmain
use modmpi
! Main BSE routine: sets up the BSE matrix and diagonalizes it.
implicit none
! local variables
integer ik,jk,a,b
integer ist,jst,i,j
integer lwork,info
real(8) t1
character(256) fname
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: epsinv(:,:,:)
complex(8), allocatable :: w(:),vl(:,:),vr(:,:)
complex(8), allocatable :: work(:)
! initialise global variables
call init0
call init1
call init2
call init3
! read density and potentials from file
call readstate
! read Fermi energy from a file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! check if system is metallic
t1=minval(abs(0.5d0-occsv(:,:)/occmax))
if (abs(t1-0.5d0).gt.0.01d0) then
  write(*,*)
  write(*,'("Error(bse): system is metallic: dielectric function will be too &
   &large")')
  write(*,'("Try using a different vkloff or reducing swidth")')
  write(*,*)
  stop
end if
! generate the BSE state index arrays
call genidxbse
! allocate global BSE arrays
if (allocated(hmlbse)) deallocate(hmlbse)
allocate(hmlbse(nmbse,nmbse))
if (allocated(evalbse)) deallocate(evalbse)
allocate(evalbse(nmbse))
! read in the RPA inverse dielectric function for q = 0
allocate(epsinv(ngrpa,ngrpa,nwrpa))
fname='EPSINV_RPA.OUT'
call getcf2pt(fname,vql(:,iq0),ngrpa,nwrpa,epsinv)
! synchronise MPI processes
call mpi_barrier(mpi_comm_world,ierror)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(bse): setting up BSE Hamiltonian matrix")')
end if
! zero the BSE Hamiltonian
hmlbse(:,:)=0.d0
! compute diagonal matrix elements
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  jk=ikmap(ivk(1,ik),ivk(2,ik),ivk(3,ik))
  do i=1,nvbse
    ist=istbse(i,ik)
    do j=1,ncbse
      jst=jstbse(j,ik)
      a=ijkbse(i,j,ik)
      hmlbse(a,a)=(evalsv(jst,jk)+scissor)-evalsv(ist,jk)
      if (bsefull) then
        b=a+nbbse
        hmlbse(b,b)=-hmlbse(a,a)
      end if
    end do
  end do
end do
deallocate(epsinv)
! add the exchange matrix elements
if (hxbse) call hmlxbse
! generate the exp(iG.r) functions for all the MBPT G-vectors
call genexpigr
! add the direct matrix elements
if (hdbse) call hmldbse
deallocate(expgmt,expgir)
! add matrices from all processes and redistribute
if (np_mpi.gt.1) then
  call mpi_allreduce(mpi_in_place,hmlbse,nmbse*nmbse,mpi_double_complex, &
   mpi_sum,mpi_comm_world,ierror)
end if
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(bse): diagonalising the BSE Hamiltonian matrix")')
end if
! diagonalize the BSE matrix
if (bsefull) then
! full non-Hermitian matrix
  allocate(w(nmbse))
  allocate(vl(1,1),vr(nmbse,nmbse))
  lwork=2*nmbse
  allocate(rwork(lwork))
  allocate(work(lwork))
  call zgeev('N','V',nmbse,hmlbse,nmbse,w,vl,1,vr,nmbse,work,lwork,rwork,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(bse): diagonalisation failed")')
    write(*,'(" ZGEEV returned INFO = ",I8)') info
    write(*,*)
    stop
  end if
  evalbse(:)=dble(w(:))
  hmlbse(:,:)=vr(:,:)
  deallocate(vl,vr,rwork,work)
else
! Hermitian block only
  allocate(rwork(3*nmbse))
  lwork=2*nmbse
  allocate(work(lwork))
  call zheev('V','U',nmbse,hmlbse,nmbse,evalbse,work,lwork,rwork,info)
  if (info.ne.0) then
    write(*,*)
    write(*,'("Error(bse): diagonalisation failed")')
    write(*,'(" ZHEEV returned INFO = ",I8)') info
    write(*,*)
    stop
  end if
  deallocate(rwork,work)
end if
! write the BSE eigenvalues to file
if (mp_mpi) then
  open(50,file='EIGVAL_BSE.OUT',action='WRITE',form='FORMATTED')
  write(50,'(I6," : nmbse")') nmbse
  if (bsefull) then
    do a=1,nmbse
      write(50,'(I6,2G18.10)') a,dble(w(a)),aimag(w(a))
    end do
    deallocate(w)
  else
    do a=1,nmbse
      write(50,'(I6,G18.10)') a,evalbse(a)
    end do
  end if
  close(50)
end if
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(bse): calculating the macroscopic dielectric tensor")')
end if
! calculate the macroscopic dielectric tensor
if (mp_mpi) call dielectric_bse
! deallocate global MBPT and BSE arrays
deallocate(istbse,jstbse,ijkbse,hmlbse,evalbse)
return
end subroutine

