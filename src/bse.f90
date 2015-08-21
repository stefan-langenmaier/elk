
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
integer ist,jst,i,j,k
integer ntop,lwork,info
real(8) h0,t1
character(256) fname
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: epsinv(:,:,:)
complex(8), allocatable :: w(:)
complex(8), allocatable :: vl(:,:),vr(:,:)
complex(8), allocatable :: work(:)
! initialise global variables
call init0
call init1
call init2
call init3
! check if the BSE extra valence or conduction states are in range
do i=1,nvxbse
  ist=istxbse(i)
  if ((ist.lt.1).or.(ist.gt.nstsv)) then
    write(*,*)
    write(*,'("Error(bse): extra valence state out of range : ",I8)') ist
    write(*,*)
    stop
  end if
end do
do j=1,ncxbse
  jst=jstxbse(j)
  if ((jst.lt.1).or.(jst.gt.nstsv)) then
    write(*,*)
    write(*,'("Error(bse): extra conduction state out of range : ",I8)') jst
    write(*,*)
    stop
  end if
end do
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
! number of valence states for transitions
nvbse=nvbse0+nvxbse
! number of conduction states for transitions
ncbse=ncbse0+ncxbse
if ((nvbse.le.0).or.(ncbse.le.0)) then
  write(*,*)
  write(*,'("Error(bse): invalid number of valence or conduction transition &
   &states : ",2I8)') nvbse,ncbse
  write(*,*)
  stop
end if
! total number of transitions
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
  if ((ntop-nvbse0+1).lt.1) then
    write(*,*)
    write(*,'("Error(bse): not enough valence states, reduce nvbse")')
    write(*,*)
    stop
  end if
  if ((ntop+ncbse0).gt.nstsv) then
    write(*,*)
    write(*,'("Error(bse): not enough conduction states, reduce ncbse or &
     &increase nempty")')
    write(*,*)
    stop
  end if
! index from BSE valence states to second-variational state numbers
  do i=1,nvbse0
    istbse(i,ik)=idx(ntop-nvbse0+i)
  end do
! index from BSE conduction states to second-variational state numbers
  do j=1,ncbse0
    jstbse(j,ik)=idx(ntop+j)
  end do
! add extra states to the list
  do i=1,nvxbse
    ist=istxbse(i)
    if (evalsv(ist,jk).gt.efermi) then
      write(*,*)
      write(*,'("Error(bse): extra valence state above Fermi energy : ",I6)') &
       ist
      write(*,'(" for k-point ",I8)') jk
      write(*,*)
      stop
    end if
    do k=1,nvbse0+i-1
      if (ist.eq.istbse(k,ik)) then
        write(*,*)
        write(*,'("Error(bse): redundant extra valence state : ",I6)') ist
        write(*,'(" for k-point ",I8)') jk
        write(*,*)
        stop
      end if
    end do
    istbse(nvbse0+i,ik)=ist
  end do
  do j=1,ncxbse
    jst=jstxbse(j)
    if (evalsv(jst,jk).lt.efermi) then
      write(*,*)
      write(*,'("Error(bse): extra conduction state below Fermi energy : ",&
       &I6)') jst
      write(*,'(" for k-point ",I8)') jk
      write(*,*)
      stop
    end if
    do k=1,ncbse0+j-1
      if (jst.eq.jstbse(k,ik)) then
        write(*,*)
        write(*,'("Error(bse): redundant extra conduction state : ",I6)') jst
        write(*,'(" for k-point ",I8)') jk
        write(*,*)
        stop
      end if
    end do
    jstbse(ncbse0+j,ik)=jst
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
! read in the RPA inverse dielectric function for q = 0
allocate(epsinv(ngrpa,ngrpa,nwrpa))
fname='EPSINV_RPA.OUT'
call getcf2pt(fname,vql(:,iq0),ngrpa,nwrpa,epsinv)
! compute the G = G' = q = 0 part of the direct kernel
h0=-2.d0/twopi**2
h0=h0*(6.d0*pi**2*(wkptnr/omega))**(1.d0/3.d0)
h0=h0*fourpi*dble(epsinv(1,1,1))
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
      hmlbse(a,a)=(evalsv(jst,jk)+scissor)-evalsv(ist,jk)+h0
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
  write(*,'("Info(bse): diagonalizing the BSE Hamiltonian matrix")')
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

