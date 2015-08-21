
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
integer ntop,lwork,info
real(8) t1
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: epsinv(:,:,:)
complex(8), allocatable :: work(:)
if (bsefull) then
  write(*,*)
  write(*,'("Error(bse): full BSE matrix not yet implimented")')
  write(*,*)
  stop
end if
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
if ((abs(t1-0.5d0).gt.0.01d0).and.mp_mpi) then
  write(*,*)
  write(*,'("Error(bse): system is metallic: dielectric function will be too &
   &large")')
  write(*,'("Try using a different vkloff or reducing swidth")')
  write(*,*)
  stop
end if
! number of transitions
nvcbse=nvbse*ncbse
! block size in BSE matrix
nbbse=nvcbse*nkptnr
! BSE matrix size
if (bsefull) then
  nmbse=2*nbbse
else
  nmbse=nbbse
end if
! allocate global BSE arrays
if (allocated(istbse)) deallocate(istbse)
allocate(istbse(nvbse,nkptnr))
if (allocated(jstbse)) deallocate(jstbse)
allocate(jstbse(ncbse,nkptnr))
if (allocated(ijkbse)) deallocate(ijkbse)
allocate(ijkbse(nvbse,ncbse,nkptnr))
if (allocated(hmlbse)) deallocate(hmlbse)
allocate(hmlbse(nmbse,nmbse))
if (allocated(evalbse)) deallocate(evalbse)
allocate(evalbse(nmbse))
allocate(idx(nstsv))
a=0
! loop over non-reduced k-points
do ik=1,nkptnr
! equivalent reduced k-point
  jk=ikmap(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! index for sorting the eigenvalues into ascending order
  call sortidx(nstsv,evalsv(:,jk),idx)
! find the topmost occupied band
  ntop=nstsv
  do ist=nstsv,1,-1
    if (evalsv(idx(ist),jk).lt.efermi) then
      ntop=ist
      exit
    end if
  end do
  if ((ntop-nvbse+1).lt.1) then
    write(*,*)
    write(*,'("Error(bse): not enough valence states, reduce nvbse")')
    write(*,*)
    stop
  end if
  if ((ntop+ncbse).gt.nstsv) then
    write(*,*)
    write(*,'("Error(bse): not enough conduction states, reduce ncbse or &
     &increase nempty")')
    write(*,*)
    stop
  end if
! index from BSE valence states to second-variational state numbers
  do i=1,nvbse
    istbse(i,ik)=idx(ntop-nvbse+i)
  end do
! index from BSE conduction states to second-variational state numbers
  do j=1,ncbse
    jstbse(j,ik)=idx(ntop+j)
  end do
! index from BSE valence-conduction pair and k-point to location in BSE matrix
  do i=1,nvbse
    do j=1,ncbse
      a=a+1
      ijkbse(i,j,ik)=a
    end do
  end do
! end loop over non-reduced k-points
end do
deallocate(idx)
! read in the RPA inverse dielectric function
allocate(epsinv(nwrpa,ngrpa,ngrpa))
call getepsinv_rpa(vql(:,iq0),epsinv)
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
call hmlxbse
! generate the exp(iG.r) functions for all the RPA G vectors
call genexpigr
! add the direct matrix elements
call hmldbse
deallocate(expgmt,expgir)
! add matrices from all processes and redistribute
if (np_mpi.gt.1) then
  call mpi_allreduce(mpi_in_place,hmlbse,nmbse*nmbse,mpi_double_complex, &
   mpi_sum,mpi_comm_world,ierror)
end if
! diagonalize the BSE matrix
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
! write the BSE eigenvalues to file
open(50,file='EIGVAL_BSE.OUT',action='WRITE',form='FORMATTED')
write(50,'(I6," : nmbse")') nmbse
do a=1,nmbse
  write(50,'(I6,G18.10)') a,evalbse(a)
end do
close(50)
! calculate the macroscopic dielectric tensor
call dielectric_bse
! deallocate global RPA and BSE arrays
deallocate(istbse,jstbse,ijkbse,hmlbse,evalbse)
return
end subroutine

