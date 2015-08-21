
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genstress
use modmain
use modphonon
use modmpi
use modstore
implicit none
! local variables
logical done(3)
integer isym,lspl,i,j
real(8) et1,s(3,3),v(3),t1
real(8) a(3,3),b(3,3),c(3,3),d(3,3)
! store original parameters
avec0(:,:)=avec(:,:)
tshift0=tshift
tshift=.false.
tforce0=tforce
tforce=.false.
done(:)=.false.
do j=1,3
  if (done(j)) cycle
  do i=1,3
    if (mp_mpi) then
      write(*,'("Info(genstress): stress matrix component : ",2I2)') i,j
    end if
    avec(:,:)=avec0(:,:)
! displace lattice vector
    avec(i,j)=avec(i,j)+deltaph
! run the ground-state calculation
    call gndstate
! subsequent calculations will read STATE.OUT
    trdstate=.true.
! store the total energy for the first displacement
    et1=engytot
! displace the lattice vector again
    avec(i,j)=avec(i,j)+deltaph
! run the ground-state calculation again
    call gndstate
! compute the stress matrix element
    a(i,j)=(engytot-et1)/deltaph
  end do
  done(j)=.true.
! copy to other stress vectors if possible
  avec(:,:)=avec0(:,:)
  call symmetry
  do isym=1,nsymcrys
    lspl=lsplsymc(isym)
    call r3mv(symlatc(:,:,lspl),avec(:,j),v)
    do i=1,3
      if (done(i)) cycle
      t1=sum(abs(avec(:,i)-v(:)))
      if (t1.lt.epslat) then
        call r3mv(symlatc(:,:,lspl),a(:,j),a(:,i))
        done(i)=.true.
      end if
    end do
  end do
end do
! restore original parameters
avec(:,:)=avec0(:,:)
tshift=tshift0
tforce=tforce0
! regenerate symmetries
call symmetry
! symmetrise the stress matrix
b(:,:)=0.d0
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  s(:,:)=dble(symlat(:,:,lspl))
  call r3mtm(symlatc(:,:,lspl),a,c)
  call r3mm(c,s,d)
  b(:,:)=b(:,:)+d(:,:)
end do
stress(:,:)=b(:,:)/dble(nsymcrys)
! compute the maximum stress magnitude over all latice vectors
stressmax=0.d0
do i=1,3
  t1=sqrt(stress(1,i)**2+stress(2,i)**2+stress(3,i)**2)
  if (t1.gt.stressmax) stressmax=t1
end do
return
end subroutine

