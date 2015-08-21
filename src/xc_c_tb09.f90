
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xc_c_tb09
use modmain
implicit none
! local variables
integer is,ias,nr,ir,itp,i
real(8), parameter :: alpha=-0.012d0, beta=1.023d0
real(8) t1
! allocatable arrays
real(8), allocatable :: grfmt(:,:,:,:),grfir(:,:)
real(8), allocatable :: rfmt(:,:,:),rfir(:)
real(8), allocatable :: rfmt1(:,:),rfmt2(:,:,:)
! external functions
real(8) rfint
external rfint
! compute the gradient of the density
allocate(grfmt(lmmaxvr,nrmtmax,natmtot,3))
allocate(grfir(ngrtot,3))
call gradrf(rhomt,rhoir,grfmt,grfir)
allocate(rfmt(lmmaxvr,nrmtmax,natmtot),rfir(ngrtot))
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt1,rfmt2,is,nr) &
!$OMP PRIVATE(i,ir,itp,t1)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt1(lmmaxvr,nrmtmax),rfmt2(lmmaxvr,nrmtmax,3))
  is=idxis(ias)
  nr=nrmt(is)
! convert muffin-tin density to spherical coordinates
  call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rhomt(:,:,ias), &
   lmmaxvr,0.d0,rfmt1,lmmaxvr)
! convert muffin-tin gradient to spherical coordinates
  do i=1,3
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
     grfmt(:,:,ias,i),lmmaxvr,0.d0,rfmt2(:,:,i),lmmaxvr)
  end do
! integrand in muffin-tin
  do ir=1,nr
    do itp=1,lmmaxvr
      t1=sqrt(rfmt2(itp,ir,1)**2+rfmt2(itp,ir,2)**2+rfmt2(itp,ir,3)**2)
      rfmt1(itp,ir)=t1/rfmt1(itp,ir)
    end do
  end do
! convert to spherical harmonics
  call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rfmt1,lmmaxvr, &
   0.d0,rfmt(:,:,ias),lmmaxvr)
  deallocate(rfmt1,rfmt2)
end do
!$OMP END DO
!$OMP END PARALLEL
! integrand in interstitial
do ir=1,ngrtot
  t1=sqrt(grfir(ir,1)**2+grfir(ir,2)**2+grfir(ir,3)**2)
  rfir(ir)=t1/rhoir(ir)
end do
! integrate over the unit cell
t1=rfint(rfmt,rfir)
! set the constant
c_tb09=alpha+beta*sqrt(abs(t1)/omega)
deallocate(grfmt,grfir,rfmt,rfir)
return
end subroutine

