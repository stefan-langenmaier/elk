
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phgvec
use modmain
use modphonon
implicit none
! local variables
integer ig,i,j
real(8) v1(3),v2(3),t1
! automatic arrays
real(8) vc(3,nphsc)
if ((task.ne.200).and.(task.ne.201)) return
if (nphsc.lt.3) return
! allocate global G-vector phonon classification array
if (allocated(igph)) deallocate(igph)
allocate(igph(ngvec))
! vectors representing each G-vector class
do i=1,nphsc
  vc(:,i)=dble(i-1)*vqc(:,iqph)
end do
! loop over G-vectors
do ig=1,ngvec
! loop over classes
  do i=1,nphsc
    v1(:)=vgc(:,ig)-vc(:,i)
    v2(:)=binv0(:,1)*v1(1)+binv0(:,2)*v1(2)+binv0(:,3)*v1(3)
    do j=1,3
      t1=v2(j)-nint(v2(j))
      if (abs(t1).gt.epslat) goto 10
    end do
! assign G-vector to class
    igph(ig)=i
    goto 20
10 continue
  end do
  write(*,*)
  write(*,'("Error(phgvec): cannot classify G-vector : ",3G18.10)') vgc(:,ig)
  write(*,*)
  stop
20 continue
end do
return
end subroutine

