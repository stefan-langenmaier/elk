
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine seceqnit(nmatp,ngp,igpig,vpl,vgpl,vgpc,apwalm,evalfv,evecfv)
use modmain
implicit none
! arguments
integer, intent(in) :: nmatp
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vgpl(3,ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer ist,jst,ias,it
real(8) ts1,ts0
real(8) t1
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: h(:),o(:)
complex(8), allocatable :: hv(:,:),ov(:,:)
! external functions
complex(8) zdotc
external zdotc
if ((iscl.ge.2).or.trdstate) then
! read in the eigenvalues/vectors from file
  call getevalfv(vpl,evalfv)
  call getevecfv(vpl,vgpl,evecfv)
else
! initialise the eigenvectors to canonical basis vectors
  evecfv(1:nmatp,:)=0.d0
  do ist=1,nstfv
    evecfv(ist,ist)=1.d0
  end do
end if
! compute Hamiltonian and overlap matrices
call timesec(ts0)
allocate(h(nmatp**2),o(nmatp**2))
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) PRIVATE(ias)
!$OMP SECTION
! Hamiltonian
h(:)=0.d0
do ias=1,natmtot
  call hmlaa(ias,ngp,apwalm,h)
  call hmlalo(ias,ngp,apwalm,h)
  call hmllolo(ias,ngp,h)
end do
call hmlistl(ngp,igpig,vgpc,h)
!$OMP SECTION
! overlap
o(:)=0.d0
do ias=1,natmtot
  call olpaa(ias,ngp,apwalm,o)
  call olpalo(ias,ngp,apwalm,o)
  call olplolo(ias,ngp,o)
end do
call olpistl(ngp,igpig,o)
!$OMP END PARALLEL SECTIONS
call timesec(ts1)
!$OMP CRITICAL
timemat=timemat+ts1-ts0
!$OMP END CRITICAL
call timesec(ts0)
allocate(hv(nmatp,nstfv),ov(nmatp,nstfv))
! start iteration loop
do it=1,nseqit
! operate with H and O on the current vectors
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ist=1,nstfv
    call zhemv('U',nmatp,zone,h,nmatp,evecfv(:,ist),1,zzero,hv(:,ist),1)
  end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ist=1,nstfv
    call zhemv('U',nmatp,zone,o,nmatp,evecfv(:,ist),1,zzero,ov(:,ist),1)
  end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(t1)
!$OMP DO
  do ist=1,nstfv
! normalise
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,ov(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      call zdscal(nmatp,t1,evecfv(:,ist),1)
      call zdscal(nmatp,t1,hv(:,ist),1)
      call zdscal(nmatp,t1,ov(:,ist),1)
    end if
! estimate the eigenvalue
    evalfv(ist)=dble(zdotc(nmatp,evecfv(:,ist),1,hv(:,ist),1))
! subtract the gradient of the Rayleigh quotient from the eigenvector
    t1=evalfv(ist)
    evecfv(1:nmatp,ist)=evecfv(1:nmatp,ist)-tauseq*(hv(1:nmatp,ist) &
     -t1*ov(1:nmatp,ist))
  end do
!$OMP END DO
!$OMP END PARALLEL
! normalise again
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ist=1,nstfv
    call zhemv('U',nmatp,zone,o,nmatp,evecfv(:,ist),1,zzero,ov(:,ist),1)
  end do
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(t1)
!$OMP DO
  do ist=1,nstfv
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,ov(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      call zdscal(nmatp,t1,evecfv(:,ist),1)
      call zdscal(nmatp,t1,ov(:,ist),1)
    end if
  end do
!$OMP END DO
!$OMP END PARALLEL
! perform Gram-Schmidt orthonormalisation
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jst,zt1,t1)
!$OMP DO ORDERED
  do ist=1,nstfv
!$OMP ORDERED
    do jst=1,ist-1
      zt1=-zdotc(nmatp,evecfv(:,jst),1,ov(:,ist),1)
      call zaxpy(nmatp,zt1,evecfv(:,jst),1,evecfv(:,ist),1)
      call zaxpy(nmatp,zt1,ov(:,jst),1,ov(:,ist),1)
    end do
!$OMP END ORDERED
! normalise
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,ov(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      call zdscal(nmatp,t1,evecfv(:,ist),1)
      call zdscal(nmatp,t1,ov(:,ist),1)
    end if
  end do
!$OMP END DO
!$OMP END PARALLEL
! end iteration loop
end do
deallocate(h,o,hv,ov)
call timesec(ts1)
!$OMP CRITICAL
timefv=timefv+ts1-ts0
!$OMP END CRITICAL
return
end subroutine

