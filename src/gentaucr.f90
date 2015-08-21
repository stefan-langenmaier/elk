
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentaucr(taumt)
use modmain
implicit none
! arguments
real(8), intent(inout) :: taumt(lmmaxvr,nrmtmax,natmtot)
! local variables
integer is,ia,ias,ist,m,ispn
integer nr,nrc,ir,irc,itp
real(8) t0,t1
! allocatable arrays
complex(8), allocatable :: wfcr(:,:,:)
complex(8), allocatable :: zfmt(:,:)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfcr,zfmt,is,ia) &
!$OMP PRIVATE(nr,nrc,ist,m,t0,t1) &
!$OMP PRIVATE(ispn,irc,ir,itp)
!$OMP DO
do ias=1,natmtot
  allocate(wfcr(lmmaxvr,nrcmtmax,2))
  allocate(zfmt(lmmaxvr,nrcmtmax))
  is=idxis(ias)
  ia=idxia(ias)
  nr=nrmt(is)
  nrc=nrcmt(is)
  do ist=1,spnst(is)
    if (spcore(ist,is)) then
      do m=-spk(ist,is),spk(ist,is)-1
        t0=evalcr(ist,ias)
! generate the core wavefunction
        call wavefcr(.true.,lradstp,is,ia,ist,m,nrcmtmax,wfcr)
        do ispn=1,2
! convert to spherical coordinates
          call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zbshtvr,lmmaxvr, &
           wfcr(:,:,ispn),lmmaxvr,zzero,zfmt,lmmaxvr)
          irc=0
          do ir=1,nr,lradstp
            irc=irc+1
            do itp=1,lmmaxvr
              t1=dble(zfmt(itp,irc))**2+aimag(zfmt(itp,irc))**2
              taumt(itp,ir,ias)=taumt(itp,ir,ias)+t0*t1
            end do
          end do
        end do
      end do
    end if
  end do
  deallocate(wfcr,zfmt)
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

