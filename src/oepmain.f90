
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepmain
use modmain
use modmpi
implicit none
! local variables
integer is,ia,ias,ik
integer ir,irc,idm
integer ld,it,n
real(8) tau,resp,t1
! allocatable arrays
real(8), allocatable :: rfmt(:,:,:)
real(8), allocatable :: rfir(:)
real(8), allocatable :: rvfmt(:,:,:,:)
real(8), allocatable :: rvfir(:,:)
real(8), allocatable :: dvxmt(:,:,:)
real(8), allocatable :: dvxir(:)
real(8), allocatable :: dbxmt(:,:,:,:)
real(8), allocatable :: dbxir(:,:)
complex(8), allocatable :: vnlcv(:,:,:,:)
complex(8), allocatable :: vnlvv(:,:,:)
! external functions
real(8) rfinp
external rfinp
if (iscl.lt.1) return
ld=lmmaxvr*lradstp
! calculate nonlocal matrix elements
allocate(vnlcv(ncrmax,natmtot,nstsv,nkpt))
allocate(vnlvv(nstsv,nstsv,nkpt))
call oepvnl(vnlcv,vnlvv)
! allocate local arrays
allocate(rfmt(lmmaxvr,nrmtmax,natmtot))
allocate(rfir(ngrtot))
allocate(dvxmt(lmmaxvr,nrcmtmax,natmtot))
allocate(dvxir(ngrtot))
if (spinpol) then
  allocate(rvfmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(rvfir(ngrtot,ndmag))
  allocate(dbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
  allocate(dbxir(ngrtot,ndmag))
end if
! set the exchange potential to zero
zvxmt(:,:,:)=0.d0
zvxir(:)=0.d0
if (spinpol) then
  zbxmt(:,:,:,:)=0.d0
  zbxir(:,:)=0.d0
end if
resp=0.d0
! initial step size
tau=tauoep(1)
!------------------------------!
!     start iteration loop     !
!------------------------------!
do it=1,maxitoep
  if (mod(it,10).eq.0) then
    write(*,'("Info(oepmain): done ",I4," iterations of ",I4)') it,maxitoep
  end if
! zero the residuals
  dvxmt(:,:,:)=0.d0
  dvxir(:)=0.d0
  if (spinpol) then
    dbxmt(:,:,:,:)=0.d0
    dbxir(:,:)=0.d0
  end if
! calculate the k-dependent residuals
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ik=1,nkpt
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    call oepresk(ik,vnlcv,vnlvv,dvxmt,dvxir,dbxmt,dbxir)
  end do
!$OMP END DO
!$OMP END PARALLEL
! add residuals from each process and redistribute
  if (np_mpi.gt.1) then
    n=lmmaxvr*nrcmtmax*natmtot
    call mpi_allreduce(mpi_in_place,dvxmt,n,mpi_double_precision,mpi_sum, &
     mpi_comm_world,ierror)
    call mpi_allreduce(mpi_in_place,dvxir,ngrtot,mpi_double_precision, &
     mpi_sum,mpi_comm_world,ierror)
    if (spinpol) then
      n=n*ndmag
      call mpi_allreduce(mpi_in_place,dbxmt,n,mpi_double_precision,mpi_sum, &
       mpi_comm_world,ierror)
      n=ngrtot*ndmag
      call mpi_allreduce(mpi_in_place,dbxir,n,mpi_double_precision,mpi_sum, &
       mpi_comm_world,ierror)
    end if
  end if
! convert muffin-tin residuals to spherical harmonics
  do is=1,nspecies
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ias,idm) SHARED(is)
!$OMP DO
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rfshtvr,lmmaxvr, &
       dvxmt(:,:,ias),lmmaxvr,0.d0,rfmt(:,:,ias),ld)
      do idm=1,ndmag
        call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rfshtvr,lmmaxvr, &
         dbxmt(:,:,ias,idm),lmmaxvr,0.d0,rvfmt(:,:,ias,idm),ld)
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do
! symmetrise the residuals
  call symrf(lradstp,rfmt,dvxir)
  if (spinpol) call symrvf(lradstp,rvfmt,dbxir)
! magnitude of residuals
  resoep=sqrt(abs(rfinp(lradstp,rfmt,rfmt,dvxir,dvxir)))
  do idm=1,ndmag
    t1=rfinp(lradstp,rvfmt(:,:,:,idm),rvfmt(:,:,:,idm),dbxir(:,idm), &
     dbxir(:,idm))
    resoep=resoep+sqrt(abs(t1))
  end do
  resoep=resoep/omega
! adjust step size
  if (it.gt.1) then
    if (resoep.gt.resp) then
      tau=tau*tauoep(2)
    else
      tau=tau*tauoep(3)
    end if
  end if
  resp=resoep
! update complex potential and field
  do is=1,nspecies
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ias,irc,ir,idm)
!$OMP DO
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
! convert residual to spherical coordinates and subtract from complex potential
        call dgemv('N',lmmaxvr,lmmaxvr,-tau,rbshtvr,lmmaxvr,rfmt(:,ir,ias),1, &
         1.d0,zvxmt(:,irc,ias),2)
        do idm=1,ndmag
          call dgemv('N',lmmaxvr,lmmaxvr,-tau,rbshtvr,lmmaxvr, &
           rvfmt(:,ir,ias,idm),1,1.d0,zbxmt(:,irc,ias,idm),2)
        end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do
  zvxir(:)=zvxir(:)-tau*dvxir(:)
  do idm=1,ndmag
    zbxir(:,idm)=zbxir(:,idm)-tau*dbxir(:,idm)
  end do
! end iteration loop
end do
! generate the real potential and field
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
! convert to real spherical harmonics
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,zvxmt(:,irc,ias),2, &
       0.d0,rfmt(:,ir,ias),1)
      do idm=1,ndmag
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr, &
         zbxmt(:,irc,ias,idm),2,0.d0,rvfmt(:,ir,ias,idm),1)
      end do
    end do
  end do
end do
! convert potential and field from a coarse to a fine radial mesh
call rfmtctof(rfmt)
do idm=1,ndmag
  call rfmtctof(rvfmt(:,:,:,idm))
end do
! add to existing (density derived) correlation potential and field
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      vxcmt(:,ir,ias)=vxcmt(:,ir,ias)+rfmt(:,ir,ias)
      do idm=1,ndmag
        bxcmt(:,ir,ias,idm)=bxcmt(:,ir,ias,idm)+rvfmt(:,ir,ias,idm)
      end do
    end do
  end do
end do
vxcir(:)=vxcir(:)+dble(zvxir(:))
do idm=1,ndmag
  bxcir(:,idm)=bxcir(:,idm)+dble(zbxir(:,idm))
end do
! symmetrise the exchange potential and field
call symrf(1,vxcmt,vxcir)
if (spinpol) call symrvf(1,bxcmt,bxcir)
deallocate(rfmt,rfir,vnlcv,vnlvv)
deallocate(dvxmt,dvxir)
if (spinpol) then
  deallocate(rvfmt,rvfir)
  deallocate(dbxmt,dbxir)
end if
return
end subroutine

