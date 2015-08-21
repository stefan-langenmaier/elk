
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotks
use modmain
use modphonon
implicit none
! local variables
integer is,ias
! allocatable arrays
complex(8), allocatable :: zfmt(:,:)
! compute the exchange-correlation potential derivative
! (at this stage the density derivative is in spherical coordinates)
call dpotxc
! convert density derivative to spherical harmonics
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zfmt,is)
!$OMP DO
do ias=1,natmtot
  allocate(zfmt(lmmaxvr,nrmtmax))
  is=idxis(ias)
  call zcopy(lmmaxvr*nrmt(is),drhomt(:,:,ias),1,zfmt,1)
  call zgemm('N','N',lmmaxvr,nrmt(is),lmmaxvr,zone,zfshtvr,lmmaxvr,zfmt, &
   lmmaxvr,zzero,drhomt(:,:,ias),lmmaxvr)
  deallocate(zfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
! generate the Coulomb potential derivative
call dpotcoul
! add to the Kohn-Sham potential derivative
do ias=1,natmtot
  is=idxis(ias)
  call zaxpy(lmmaxvr*nrmt(is),zone,dvclmt(:,:,ias),1,dvsmt(:,:,ias),1)
end do
call zaxpy(ngtot,zone,dvclir,1,dvsir,1)
! remove the gradient part of the potential derivative for displaced muffin-tin
call zaxpy(lmmaxvr*nrmt(isph),zone,gvsmt,1,dvsmt(:,:,iasph),1)
return
end subroutine

