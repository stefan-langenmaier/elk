
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getcf2pt(fname,vpl,ng,m,cf)
use modmain
! arguments
character(256), intent(in) :: fname
real(8), intent(in) :: vpl(3)
integer, intent(in) :: ng
integer, intent(in) :: m
complex(8), intent(out) :: cf(ng,ng,m)
! local variables
integer isym,iq,i
integer ig,jg,igm,jgm
integer lspl,ilspl
integer recl,ng_,m_
real(8) vql_(3),si(3,3)
real(8) vgql(3),v(3),t1
! allocatable arrays
integer, allocatable :: map(:)
real(8), allocatable :: vgpl(:,:)
complex(8), allocatable :: cf_(:,:,:)
complex(8), allocatable :: zv(:)
! find the equivalent reduced q-point
call findqpt(vpl,isym,iq)
! find the record length
inquire(iolength=recl) vql(:,1),ng,m,cf
!$OMP CRITICAL
open(100,file=trim(fname),action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
read(100,rec=iq) vql_,ng_,m_,cf
close(100)
!$OMP END CRITICAL
t1=abs(vql(1,iq)-vql_(1))+abs(vql(2,iq)-vql_(2))+abs(vql(3,iq)-vql_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getcf2pt): differing vectors for q-point ",I8)') iq
  write(*,'(" current : ",3G18.10)') vql(:,iq)
  write(*,'(" file    : ",3G18.10)') vql_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
if (ng.ne.ng_) then
  write(*,*)
  write(*,'("Error(getcf2pt): differing ng for q-point ",I8)') iq
  write(*,'(" current : ",I8)') ng
  write(*,'(" file    : ",I8)') ng_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
if (m.ne.m_) then
  write(*,*)
  write(*,'("Error(getcf2pt): differing m for q-point ",I8)') iq
  write(*,'(" current : ",I8)') m
  write(*,'(" file    : ",I8)') m_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
! if p = q then return
t1=abs(vpl(1)-vql(1,iq))+abs(vpl(2)-vql(2,iq))+abs(vpl(3)-vql(3,iq))
if (t1.lt.epslat) return
! allocate local arrays
allocate(map(ng))
allocate(vgpl(3,ng))
allocate(cf_(ng,ng,m))
allocate(zv(ng))
! perform translation operation and store in temporary array
if (tvzsymc(isym)) then
! translation vector is zero
  cf_(:,:,:)=cf(:,:,:)
else
! non-zero translation vector gives a phase factor
  do ig=1,ng
    t1=twopi*dot_product(dble(ivg(:,ig)),vtlsymc(:,isym))
    zv(ig)=cmplx(cos(t1),sin(t1),8)
  end do
  do ig=1,ng
    do jg=1,ng
      cf_(ig,jg,:)=zv(ig)*conjg(zv(jg))*cf(ig,jg,:)
    end do
  end do
end if
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! the inverse of the spatial symmetry rotates q into p
ilspl=isymlat(lspl)
si(:,:)=dble(symlat(:,:,ilspl))
! find the map from {G+q} to {G+p}
map(:)=0
do ig=1,ng
  vgpl(:,ig)=dble(ivg(:,ig))+vpl(:)
end do
i=1
do ig=1,ng
  vgql(:)=dble(ivg(:,ig))+vql(:,iq)
  call r3mtv(si,vgql,v)
  do jg=i,ng
    t1=abs(v(1)-vgpl(1,jg))+abs(v(2)-vgpl(2,jg))+abs(v(3)-vgpl(3,jg))
    if (t1.lt.epslat) then
      map(ig)=jg
      if (jg.eq.i) i=i+1
      exit
    end if
  end do
end do
! rotate epsilon inverse
do ig=1,ng
  igm=map(ig)
  do jg=1,ng
    jgm=map(jg)
    if ((igm.eq.0).or.(jgm.eq.0)) then
      cf(ig,jg,:)=0.d0
    else
      cf(ig,jg,:)=cf_(igm,jgm,:)
    end if
  end do
end do
deallocate(map,vgpl,cf_,zv)
return
end subroutine

