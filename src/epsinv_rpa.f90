
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine epsinv_rpa
use modmain
use modmpi
implicit none
! local variables
integer ik,iq,ig,iw
integer ist,jst,n
integer info1,info2,recl
real(8) vgqc(3),eij
real(8) t1,t2,t3
! allocatable arrays
integer, allocatable :: ipiv(:)
complex(8), allocatable :: epsinv(:,:,:)
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: a(:,:)
complex(8), allocatable :: work(:)
! initialise global variables
call init0
call init1
call init2
call init3
allocate(epsinv(nwrpa,ngrpa,ngrpa))
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
if (mp_mpi) then
! determine the record length for EPSINV_RPA.OUT
  inquire(iolength=recl) vql(:,1),nwrpa,ngrpa,epsinv
! open EPSINV_RPA.OUT
  open(50,file='EPSINV_RPA.OUT',action='WRITE',form='UNFORMATTED', &
   access='DIRECT',status='REPLACE',recl=recl)
end if
! loop over q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(epsinv_rpa): ",I6," of ",I6," q-points")') iq,nqpt
  epsinv(:,:,:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! compute response function chi0 and store in array epsinv
    call genchi0(iq,ik,epsinv)
  end do
!$OMP END DO
!$OMP END PARALLEL
! add chi0 from each process and redistribute
  if (np_mpi.gt.1) then
    n=nwrpa*ngrpa*ngrpa
    call mpi_allreduce(mpi_in_place,epsinv,n,mpi_double_complex,mpi_sum, &
     mpi_comm_world,ierror)
  end if
!-----------------------------------!
!     compute epsilon from chi0     !
!-----------------------------------!
  do ig=1,ngrpa
    vgqc(:)=vgc(:,ig)+vqc(:,iq)
    t1=vgqc(1)**2+vgqc(2)**2+vgqc(3)**2+lmda2rpa
    t1=-fourpi/t1
    epsinv(:,ig,:)=t1*epsinv(:,ig,:)
    epsinv(:,ig,ig)=epsinv(:,ig,ig)+1.d0
  end do
!----------------------------------------!
!     compute G = G' = q = 0 element     !
!----------------------------------------!
  if (iq.eq.iq0) then
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(pmat,t1,t2,t3) &
!$OMP PRIVATE(ist,jst,eij)
!$OMP DO
    do ik=1,nkpt
! distribute among MPI processes
      if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
      allocate(pmat(3,nstsv,nstsv))
! read the momentum matrix elements from file
      call getpmat(vkl(:,ik),pmat)
! prefactor: the 1/3 is because we want Tr(epsinv)/3
      t1=(8.d0/3.d0)*pi*wkpt(ik)/omega
      do ist=1,nstsv
        do jst=1,nstsv
          t2=occsv(ist,ik)*(1.d0-occsv(jst,ik)/occmax)
          if (t2.gt.1.d-8) then
            eij=evalsv(ist,ik)-(evalsv(jst,ik)+scissor)
            if (abs(eij).gt.1.d-8) then
! sum |<i|p|j>|^2
              t3=sum(dble(pmat(:,ist,jst))**2+aimag(pmat(:,ist,jst))**2)
              t3=t1*t2*t3
!$OMP CRITICAL
              epsinv(:,1,1)=epsinv(:,1,1)-t3/(eij*(eij**2-wrpa(:)**2))
!$OMP END CRITICAL
            end if
          end if
        end do
      end do
      deallocate(pmat)
    end do
!$OMP END DO
!$OMP END PARALLEL
! add epsinv from each process and redistribute
    if (np_mpi.gt.1) then
      call mpi_allreduce(mpi_in_place,epsinv(:,1,1),nwrpa,mpi_double_complex, &
       mpi_sum,mpi_comm_world,ierror)
    end if
  end if
!-------------------------------------!
!     invert epsilon over G-space     !
!-------------------------------------!
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ipiv,a,work,info1,info2)
!$OMP DO
  do iw=1,nwrpa
    allocate(ipiv(ngrpa))
    allocate(a(ngrpa,ngrpa))
    allocate(work(ngrpa))
    a(:,:)=epsinv(iw,:,:)
    call zgetrf(ngrpa,ngrpa,a,ngrpa,ipiv,info1)
    call zgetri(ngrpa,a,ngrpa,ipiv,work,ngrpa,info2)
    if ((info1.ne.0).or.(info2.ne.0)) then
      write(*,*)
      write(*,'("Error(epsinv_rpa): unable to invert epsilon")')
      write(*,'(" for q-point ",I6)') iq
      write(*,'(" and RPA frequency ",I6)') iw
      write(*,*)
      stop
    end if
    epsinv(iw,:,:)=a(:,:)
    deallocate(ipiv,a,work)
  end do
!$OMP END DO
!$OMP END PARALLEL
! write inverse RPA epsilon to EPSINV_RPA.OUT
  if (mp_mpi) write(50,rec=iq) vql(:,iq),nwrpa,ngrpa,epsinv
! end loop over q-points
end do
if (mp_mpi) close(50)
deallocate(epsinv)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(epsinv_rpa):")')
  write(*,'(" inverse RPA dielectric function, eps^(-1)(G,G'',q,w), written to &
   &EPSINV_RPA.OUT")')
  write(*,*)
end if
return
end subroutine

