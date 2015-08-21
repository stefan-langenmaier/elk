
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentaucr(taumt)
use modmain
implicit none
! arguments
real(8), intent(inout) :: taumt(lmmaxvr,nrmtmax,natmtot,nspinor)
! local variables
integer is,ia,ias,nr,i
integer ist,m,ispn,jspn
! allocatable arrays
complex(8), allocatable :: wfcr(:,:,:)
complex(8), allocatable :: gzfmt(:,:,:),zfmt(:,:)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfcr,gzfmt,zfmt) &
!$OMP PRIVATE(is,ia,nr,ist,m) &
!$OMP PRIVATE(ispn,jspn,i)
!$OMP DO
do ias=1,natmtot
  allocate(wfcr(lmmaxvr,nrmtmax,2))
  allocate(gzfmt(lmmaxvr,nrmtmax,3),zfmt(lmmaxvr,nrmtmax))
  is=idxis(ias)
  ia=idxia(ias)
  nr=nrmt(is)
  do ist=1,spnst(is)
    if (spcore(ist,is)) then
      do m=-spk(ist,is),spk(ist,is)-1
! generate the core wavefunction in spherical harmonics
        call wavefcr(.true.,1,is,ia,ist,m,nrmtmax,wfcr)
        do ispn=1,2
          if (spinpol) then
            jspn=ispn
          else
            jspn=1
          end if
! compute the gradient
          call gradzfmt(lmaxvr,nr,spr(:,is),lmmaxvr,nrmtmax,wfcr(:,:,ispn), &
           gzfmt)
          do i=1,3
! convert gradient to spherical coordinates
            call zgemm('N','N',lmmaxvr,nr,lmmaxvr,zone,zbshtvr,lmmaxvr, &
             gzfmt(:,:,i),lmmaxvr,zzero,zfmt,lmmaxvr)
            taumt(:,1:nr,ias,jspn)=taumt(:,1:nr,ias,jspn) &
             +0.5d0*(dble(zfmt(:,1:nr))**2+aimag(zfmt(:,1:nr))**2)
          end do
        end do
      end do
    end if
  end do
  deallocate(wfcr,gzfmt,zfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

