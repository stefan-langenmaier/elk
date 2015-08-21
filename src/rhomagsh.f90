
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
integer ld,is,ia,ias
integer ir,irc,idm
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
ld=lmmaxvr*lradstp
do is=1,nspecies
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,ias,irc,ir,idm)
!$OMP DO
  do ia=1,natoms(is)
    allocate(rfmt(lmmaxvr,nrcmtmax))
    ias=idxas(ia,is)
! convert the density to spherical harmonics
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      rfmt(:,irc)=rhomt(:,ir,ias)
    end do
    call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rfshtvr,lmmaxvr,rfmt, &
     lmmaxvr,0.d0,rhomt(:,:,ias),ld)
! convert magnetisation to spherical harmonics
    if (spinpol) then
      do idm=1,ndmag
        irc=0
        do ir=1,nrmt(is),lradstp
          irc=irc+1
          rfmt(:,irc)=magmt(:,ir,ias,idm)
        end do
        call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rfshtvr,lmmaxvr, &
         rfmt,lmmaxvr,0.d0,magmt(:,:,ias,idm),ld)
      end do
    end if
    deallocate(rfmt)
  end do
!$OMP END DO
!$OMP END PARALLEL
end do
return
end subroutine
!EOC

