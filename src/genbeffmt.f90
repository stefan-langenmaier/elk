
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine genbeffmt
use modmain
implicit none
! local variables
integer ld,is,ia,ias
integer nrc,irc,ir,i
real(8) cb,cso,rm,t1
! allocatable arrays
real(8), allocatable :: vr(:),drv(:),cf(:,:)
if (.not.spinpol) return
! coupling constant of the external field (g_e/4c)
cb=gfacte/(4.d0*solsc)
! coefficient of spin-orbit coupling
cso=1.d0/(4.d0*solsc**2)
ld=lmmaxvr*lradstp
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(is,ia,nrc,i,t1) &
!$OMP PRIVATE(vr,drv,cf,irc,ir,rm)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias); ia=idxia(ias)
  nrc=nrcmt(is)
! exchange-correlation magnetic field in spherical coordinates
  do i=1,ndmag
    call dgemm('N','N',lmmaxvr,nrc,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
     bxcmt(:,:,ias,i),ld,0.d0,beffmt(:,:,ias,i),lmmaxvr)
  end do
! add the external magnetic field
  t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
  beffmt(:,1:nrc,ias,ndmag)=beffmt(:,1:nrc,ias,ndmag)+t1
  if (ncmag) then
    do i=1,2
      t1=cb*(bfcmt(i,ia,is)+bfieldc(i))
      beffmt(:,1:nrc,ias,i)=beffmt(:,1:nrc,ias,i)+t1
    end do
  end if
! spin-orbit coupling radial function
  if (spinorb) then
    allocate(vr(nrmtmax),drv(nrmtmax),cf(4,nrmtmax))
! radial derivative of the spherical part of the potential
    vr(1:nrmt(is))=veffmt(1,1:nrmt(is),ias)*y00
    call fderiv(1,nrmt(is),spr(:,is),vr,drv,cf)
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      rm=1.d0-2.d0*cso*vr(ir)
      socfr(irc,ias)=cso*drv(ir)/(spr(ir,is)*rm**2)
    end do
    deallocate(vr,drv,cf)
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

