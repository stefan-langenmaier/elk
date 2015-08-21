
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
integer is,ia,it
integer ist,jst
real(8) ts1,ts0
real(8) t1
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: h(:,:)
complex(8), allocatable :: o(:,:)
! external functions
complex(8) zdotc
external zdotc
call timesec(ts0)
allocate(h(nmatmax,nstfv))
allocate(o(nmatmax,nstfv))
if ((iscl.ge.2).or.(task.eq.1).or.(task.eq.3)) then
! read in the eigenvalues/vectors from file
  call getevalfv(vpl,evalfv)
  call getevecfv(vpl,vgpl,evecfv)
else
! initialise the eigenvectors to canonical basis vectors
  evecfv(:,:)=0.d0
  do ist=1,nstfv
    evecfv(ist,ist)=1.d0
  end do
end if
! start iteration loop
do it=1,nseqit
! operate with H and O on the current vectors
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) PRIVATE(is,ia)
!$OMP SECTION
  h(:,:)=0.d0
  do is=1,nspecies
    do ia=1,natoms(is)
      call hmlaa(.true.,is,ia,ngp,apwalm,evecfv,h)
      call hmlalo(.true.,is,ia,ngp,apwalm,evecfv,h)
      call hmllolo(.true.,is,ia,ngp,evecfv,h)
    end do
  end do
  call hmlistl(.true.,ngp,igpig,vgpc,evecfv,h)
!$OMP SECTION
  o(:,:)=0.d0
  do is=1,nspecies
    do ia=1,natoms(is)
      call olpaa(.true.,is,ia,ngp,apwalm,evecfv,o)
      call olpalo(.true.,is,ia,ngp,apwalm,evecfv,o)
      call olplolo(.true.,is,ia,ngp,evecfv,o)
    end do
  end do
  call olpistl(.true.,ngp,igpig,evecfv,o)
!$OMP END PARALLEL SECTIONS
  do ist=1,nstfv
! normalise
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,o(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      evecfv(1:nmatp,ist)=t1*evecfv(1:nmatp,ist)
      h(:,ist)=t1*h(:,ist)
      o(:,ist)=t1*o(:,ist)
    end if
! estimate the eigenvalue
    evalfv(ist)=dble(zdotc(nmatp,evecfv(:,ist),1,h(:,ist),1))
! subtract the gradient of the Rayleigh quotient from the eigenvector
    t1=evalfv(ist)
    evecfv(1:nmatp,ist)=evecfv(1:nmatp,ist)-tauseq*(h(1:nmatp,ist) &
     -t1*o(1:nmatp,ist))
  end do
! normalise again
  o(:,:)=0.d0
  do is=1,nspecies
    do ia=1,natoms(is)
      call olpaa(.true.,is,ia,ngp,apwalm,evecfv,o)
      call olpalo(.true.,is,ia,ngp,apwalm,evecfv,o)
      call olplolo(.true.,is,ia,ngp,evecfv,o)
    end do
  end do
  call olpistl(.true.,ngp,igpig,evecfv,o)
  do ist=1,nstfv
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,o(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      evecfv(1:nmatp,ist)=t1*evecfv(1:nmatp,ist)
      o(:,ist)=t1*o(:,ist)
    end if
! end loop over states
  end do
! perform Gram-Schmidt orthonormalisation
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(jst,zt1,t1)
!$OMP DO ORDERED
  do ist=1,nstfv
!$OMP ORDERED
    do jst=1,ist-1
      zt1=-zdotc(nmatp,evecfv(:,jst),1,o(:,ist),1)
      call zaxpy(nmatp,zt1,evecfv(:,jst),1,evecfv(:,ist),1)
      call zaxpy(nmatp,zt1,o(:,jst),1,o(:,ist),1)
    end do
!$OMP END ORDERED
! normalise
    t1=dble(zdotc(nmatp,evecfv(:,ist),1,o(:,ist),1))
    if (t1.gt.0.d0) then
      t1=1.d0/sqrt(t1)
      evecfv(1:nmatp,ist)=t1*evecfv(1:nmatp,ist)
      o(:,ist)=t1*o(:,ist)
    end if
  end do
!$OMP END DO
!$OMP END PARALLEL
! end iteration loop
end do
deallocate(h,o)
call timesec(ts1)
!$OMP CRITICAL
timefv=timefv+ts1-ts0
!$OMP END CRITICAL
return
end subroutine

