
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzvclmt(nr,ld1,r,ld2,zrhomt,zvclmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr(nspecies)
integer, intent(in) :: ld1
real(8), intent(in) :: r(ld1,nspecies)
integer, intent(in) :: ld2
complex(8), intent(in) :: zrhomt(lmmaxvr,ld2,natmtot)
complex(8), intent(out) :: zvclmt(lmmaxvr,ld2,natmtot)
! local variables
integer is,ia,ias
do is=1,nspecies
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ias)
!$OMP DO
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    call zpotclmt(lmaxvr,nr(is),r(:,is),lmmaxvr,zrhomt(:,:,ias),zvclmt(:,:,ias))
  end do
!$OMP END DO
!$OMP END PARALLEL
end do
return
end subroutine

