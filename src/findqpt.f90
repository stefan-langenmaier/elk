
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findqpt(vpl,isym,iq)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(out) :: isym
integer, intent(out) :: iq
! local variables
integer lspl,iv(3)
real(8) v1(3),v2(3),t1
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
! multiply transpose of symmetry matrix with vpl
  v1(:)=symlat(1,:,lspl)*vpl(1) &
       +symlat(2,:,lspl)*vpl(2) &
       +symlat(3,:,lspl)*vpl(3)
! map vector components to [0,1) interval
  call r3frac(epslat,v1,iv)
! search q-points for this vector
  do iq=1,nqpt
    v2(:)=vql(:,iq)
    call r3frac(epslat,v2,iv)
    t1=abs(v1(1)-v2(1))+abs(v1(2)-v2(2))+abs(v1(3)-v2(3))
    if (t1.lt.epslat) return
  end do
end do
write(*,*)
write(*,'("Error(findqpt): equivalent q-point not in set")')
write(*,'(" Requested q-point : ",3G18.10)') vpl
write(*,*)
stop
end subroutine

