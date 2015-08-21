
! Copyright (C) 2010 Alexey I. Baranov.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genhvec
use modmain
implicit none
! local variables
logical lsym(48)
integer nh(3),ihv(3,2),nht
integer nsym,isym,sym(3,3,48)
integer ih,jh,i1,i2,i3,iv(3),k
real(8) hm2,v(3),t1
! allocatable arrays
integer, allocatable :: idx(:)
integer, allocatable :: iar(:)
real(8), allocatable :: h2(:)
! number of H-vectors in each dimension
nh(:)=int(hmax*sqrt(avec(1,:)**2+avec(2,:)**2+avec(3,:)**2)/pi)+1
! total number of H-vectors in box
nht=nh(1)*nh(2)*nh(3)
allocate(h2(nht))
! allocate global arrays
if (allocated(ivh)) deallocate(ivh)
allocate(ivh(3,nht))
if (allocated(mulh)) deallocate(mulh)
allocate(mulh(nht))
! H-vector intervals
ihv(:,1)=nh(:)/2-nh(:)+1
ihv(:,2)=nh(:)/2
! find the subgroup of symmorphic, non-magnetic symmetries
lsym(:)=.false.
do isym=1,nsymcrys
  if (tvzsymc(isym).and.(lspnsymc(isym).eq.1)) lsym(lsplsymc(isym))=.true.
end do
nsym=0
do isym=1,nsymlat
  if (lsym(isym)) then
    nsym=nsym+1
    sym(:,:,nsym)=symlat(:,:,isym)
  end if
end do
! use only the identity element if no reduction is required
if (.not.reduceh) nsym=1
! generate the symmetry inequivalent H-vectors and multiplicities
hm2=hmax**2
ih=0
do i1=ihv(1,1),ihv(1,2)
  do i2=ihv(2,1),ihv(2,2)
    do i3=ihv(3,1),ihv(3,2)
      v(:)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
      t1=v(1)**2+v(2)**2+v(3)**2
      if (t1.lt.hm2) then
        do jh=1,ih
          if (abs(h2(jh)-t1).lt.epslat) then
            do isym=1,nsym
              call i3mtv(sym(:,:,isym),ivh(:,jh),iv)
              k=abs(i1-iv(1))+abs(i2-iv(2))+abs(i3-iv(3))
              if (k.eq.0) then
                mulh(jh)=mulh(jh)+1
                goto 10
              end if
            end do
          end if
        end do
        ih=ih+1
        h2(ih)=t1
        ivh(1,ih)=i1; ivh(2,ih)=i2; ivh(3,ih)=i3
        mulh(ih)=1
      end if
10 continue
    end do
  end do
end do
nhvec=ih
! sort the H-vectors by length
allocate(idx(nhvec),iar(nhvec))
call sortidx(nhvec,h2,idx)
do k=1,3
  iar(1:nhvec)=ivh(k,1:nhvec)
  do ih=1,nhvec
    ivh(k,ih)=iar(idx(ih))
  end do
end do
iar(1:nhvec)=mulh(1:nhvec)
do ih=1,nhvec
  mulh(ih)=iar(idx(ih))
end do
deallocate(idx,iar,h2)
return
end subroutine

