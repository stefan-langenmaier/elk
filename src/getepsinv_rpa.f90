
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getepsinv_rpa(vpl,epsinv)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: epsinv(nwrpa,ngrpa,ngrpa)
! local variables
integer isym,iq,i
integer ig,jg,igm,jgm
integer lspl,ilspl
integer recl,nwrpa_,ngrpa_
real(8) vql_(3),si(3,3)
real(8) vgql(3),v(3),t1
! allocatable arrays
integer, allocatable :: map(:)
real(8), allocatable :: vgpl(:,:)
complex(8), allocatable :: epsinvt(:,:,:)
complex(8), allocatable :: zv(:)
! find the equivalent reduced q-point
call findqpt(vpl,isym,iq)
! find the record length
inquire(iolength=recl) vql(:,1),nwrpa,ngrpa,epsinv
!$OMP CRITICAL
open(100,file='EPSINV_RPA.OUT',action='READ',form='UNFORMATTED', &
 access='DIRECT',recl=recl)
read(100,rec=iq) vql_,nwrpa_,ngrpa_,epsinv
close(100)
!$OMP END CRITICAL
t1=abs(vql(1,iq)-vql_(1))+abs(vql(2,iq)-vql_(2))+abs(vql(3,iq)-vql_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getepsinv_rpa): differing vectors for q-point ",I8)') iq
  write(*,'(" current        : ",3G18.10)') vql(:,iq)
  write(*,'(" EPSINV_RPA.OUT : ",3G18.10)') vql_
  write(*,*)
  stop
end if
if (nwrpa.ne.nwrpa_) then
  write(*,*)
  write(*,'("Error(getepsinv_rpa): differing nwrpa for q-point ",I8)') iq
  write(*,'(" current        : ",I8)') nwrpa
  write(*,'(" EPSINV_RPA.OUT : ",I8)') nwrpa_
  write(*,*)
  stop
end if
if (ngrpa.ne.ngrpa_) then
  write(*,*)
  write(*,'("Error(getepsinv_rpa): differing ngrpa for q-point ",I8)') iq
  write(*,'(" current        : ",I8)') ngrpa
  write(*,'(" EPSINV_RPA.OUT : ",I8)') ngrpa_
  write(*,*)
  stop
end if
! if p = q then return
t1=abs(vpl(1)-vql(1,iq))+abs(vpl(2)-vql(2,iq))+abs(vpl(3)-vql(3,iq))
if (abs(t1).lt.epslat) return
! allocate local arrays
allocate(map(ngrpa))
allocate(vgpl(3,ngrpa))
allocate(epsinvt(nwrpa,ngrpa,ngrpa))
allocate(zv(ngrpa))
! perform translation operation and store in temporary array
if (tvzsymc(isym)) then
! translation vector is zero
  epsinvt(:,:,:)=epsinv(:,:,:)
else
! non-zero translation vector gives a phase factor
  do ig=1,ngrpa
    t1=twopi*dot_product(dble(ivg(:,ig)),vtlsymc(:,isym))
    zv(ig)=cmplx(cos(t1),sin(t1),8)
  end do
  do ig=1,ngrpa
    do jg=1,ngrpa
      epsinvt(:,ig,jg)=zv(ig)*conjg(zv(jg))*epsinv(:,ig,jg)
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
do ig=1,ngrpa
  vgpl(:,ig)=dble(ivg(:,ig))+vpl(:)
end do
i=1
do ig=1,ngrpa
  vgql(:)=dble(ivg(:,ig))+vql(:,iq)
  call r3mtv(si,vgql,v)
  do jg=i,ngrpa
    t1=abs(v(1)-vgpl(1,jg))+abs(v(2)-vgpl(2,jg))+abs(v(3)-vgpl(3,jg))
    if (t1.lt.epslat) then
      map(ig)=jg
      if (jg.eq.i) i=i+1
      goto 10
    end if
  end do
10 continue
end do
! rotate epsilon inverse
do ig=1,ngrpa
  igm=map(ig)
  do jg=1,ngrpa
    jgm=map(jg)
    if ((igm.eq.0).or.(jgm.eq.0)) then
      epsinv(:,ig,jg)=0.d0
    else
      epsinv(:,ig,jg)=epsinvt(:,igm,jgm)
    end if
  end do
end do
deallocate(map,vgpl,epsinvt,zv)
return
end subroutine

