
! Copyright (C) 2013 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zfmtctof(zfmt)
use modmain
implicit none
! arguments
real(8), intent(inout) :: zfmt(2,lmmaxvr,nrmtmax,natmtot)
! local variables
integer ld,is,ias,lm,i
ld=2*lmmaxvr*lradstp
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is,lm,i)
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
! interpolate with a clamped spline
  do lm=1,lmmaxvr
    do i=1,2
      call rfinterp(nrcmt(is),rcmt(:,is),ld,zfmt(i,lm,1,ias),nrmt(is), &
       spr(:,is),2*lmmaxvr,zfmt(i,lm,1,ias))
    end do
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
return
end subroutine

