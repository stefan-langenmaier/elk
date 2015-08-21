
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rhomagsh
! !INTERFACE:
subroutine rhomagsh
! !USES:
use modmain
! !DESCRIPTION:
!   Converts the muffin-tin density and magnetisation from spherical coordinates
!   to a spherical harmonic expansion. See {\tt rhomagk}.
!
! !REVISION HISTORY:
!   Created January 2009 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ld,idm,is,ias
integer nr,nri,iro,ir
integer nrci,nrco,irco,irc
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
ld=lmmaxvr*lradstp
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nr,nri,iro) &
!$OMP PRIVATE(nrci,nrco,irco,irc,ir,idm)
!$OMP DO
do ias=1,natmtot
  allocate(rfmt(lmmaxvr,nrcmtmax))
  is=idxis(ias)
! number of inner/outer and fine/coarse radial mesh points
  nr=nrmt(is)
  nri=nrmtinr(is)
  iro=nri+1
  nrci=nrcmtinr(is)
  nrco=nrcmt(is)-nrci
  irco=nrci+1
! convert the density to spherical harmonics
  irc=0
  do ir=1,nr,lradstp
    irc=irc+1
    rfmt(:,irc)=rhomt(:,ir,ias)
  end do
! inner part of muffin-tin
  call dgemm('N','N',lmmaxinr,nrci,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rfmt,lmmaxvr, &
   0.d0,rhomt(:,:,ias),ld)
! outer part of muffin-tin
  call dgemm('N','N',lmmaxvr,nrco,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rfmt(:,irco), &
   lmmaxvr,0.d0,rhomt(:,iro,ias),ld)
! zero the components with lmax > lmaxinr for the inner part
  do ir=1,nri,lradstp
    rhomt(lmmaxinr+1:,ir,ias)=0.d0
  end do
! convert magnetisation to spherical harmonics
  if (spinpol) then
    do idm=1,ndmag
      irc=0
      do ir=1,nr,lradstp
        irc=irc+1
        rfmt(:,irc)=magmt(:,ir,ias,idm)
      end do
! inner part of muffin-tin
      call dgemm('N','N',lmmaxinr,nrci,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rfmt, &
       lmmaxvr,0.d0,magmt(:,:,ias,idm),ld)
! outer part of muffin-tin
      call dgemm('N','N',lmmaxvr,nrco,lmmaxvr,1.d0,rfshtvr,lmmaxvr, &
       rfmt(:,irco),lmmaxvr,0.d0,magmt(:,iro,ias,idm),ld)
! zero the components with lmax > lmaxinr for the inner part
      do ir=1,nri,lradstp
        magmt(lmmaxinr+1:,ir,ias,idm)=0.d0
      end do
    end do
  end if
  deallocate(rfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine
!EOC

