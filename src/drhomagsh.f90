
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine drhomagsh
use modmain
use modphonon
implicit none
! local variables
integer is,ia,ias
integer nrc,idm
! allocatable arrays
complex(8), allocatable :: zfmt(:,:)
do is=1,nspecies
  nrc=nrcmt(is)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,ias,idm)
!$OMP DO
  do ia=1,natoms(is)
    allocate(zfmt(lmmaxvr,nrcmtmax))
    ias=idxas(ia,is)
! convert the density derivative to spherical harmonics
    zfmt(:,1:nrc)=drhomt(:,1:nrc,ias)
    call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr,zfmt,lmmaxvr, &
     zzero,drhomt(:,:,ias),lmmaxvr)
! convert the magnetisation derivative to spherical harmonics
    if (spinpol) then
      do idm=1,ndmag
        zfmt(:,1:nrc)=dmagmt(:,1:nrc,ias,idm)
        call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr,zfmt, &
         lmmaxvr,zzero,dmagmt(:,:,ias,idm),lmmaxvr)
      end do
    end if
    deallocate(zfmt)
  end do
!$OMP END DO
!$OMP END PARALLEL
end do
return
end subroutine

