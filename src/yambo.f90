
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine yambo
use modmain
implicit none
! local variables
integer nk,ik,jk,ist
integer ng_vec,ig,igk
integer isym,is,ia,xc
real(8) v(3),cv,t
character(3) symb
! allocatable arrays
logical, allocatable :: usek(:)
logical, allocatable :: useg(:)
integer, allocatable :: idx(:)
integer, allocatable :: wf_igk(:)
real(8), allocatable :: vgl(:,:)
complex(8), allocatable :: wfpw(:,:,:)
if (reducek.ne.2) then
  write(*,*)
  write(*,'("Warning(yambo): reducek should be 2 for writing Yambo input")')
  write(*,*)
end if
! initialise global variables
call init0
call init1
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! allocate local arrays
allocate(usek(nkpt))
allocate(useg(ngvec))
allocate(idx(ngvec))
allocate(wf_igk(ngkmax))
allocate(wfpw(ngkmax,nspinor,nstsv))
! remove k-points related by time-reversal symmetry
usek(:)=.true.
if (.not.spinpol) then
  do ik=1,nkpt
    if (usek(ik)) then
      v(:)=-vkl(:,ik)
      call findkpt(v,isym,jk)
      if (ik.ne.jk) usek(jk)=.false.
    end if
  end do
end if
nk=0
do ik=1,nkpt
  if (usek(ik)) nk=nk+1
end do
! find the set of G-vectors used for the wavefunctions
useg(:)=.false.
do ik=1,nkpt
  if (usek(ik)) then
    do igk=1,ngk(1,ik)
      useg(igkig(igk,1,ik))=.true.
    end do
  end if
end do
! write the G-vectors to g_vec (this is ordered since ivg is ordered)
ng_vec=0
do ig=1,ngvec
  if (useg(ig)) then
    ng_vec=ng_vec+1
    idx(ig)=ng_vec
  end if
end do
allocate(vgl(3,ng_vec))
do ig=1,ngvec
  if (useg(ig)) vgl(:,idx(ig))=dble(ivg(:,ig))
end do
! exchange-correlation type following the ABINIT convention
select case(xctype(1))
case(3)
  xc=7
case(4)
  xc=6
case(20)
  xc=11
case(21)
  xc=14
case(26)
  xc=23
case default
  write(*,*)
  write(*,'("Error(yambo): unsupported xctype : ",I6)') xctype(1)
  write(*,*)
  stop
end select
! write all necessary variables to YAMBO.OUT
open(50,file='YAMBO.OUT',action='WRITE',form='UNFORMATTED')
write(50) version
write(50) avec
write(50) nspecies
do is=1,nspecies
  symb=trim(spsymb(is))
  write(50) symb
  write(50) spzn(is)
! determine the valence charge of the atom
  cv=0.d0
  do ist=1,spnst(is)
    if (.not.spcore(ist,is)) cv=cv+spocc(ist,is)
  end do
  write(50) cv
  write(50) natoms(is)
  do ia=1,natoms(is)
    write(50) atposl(:,ia,is)
  end do
end do
write(50) chgexs
write(50) xc
! non-zero temperature only for Fermi-Dirac smearing
if (stype.eq.3) then
  t=swidth/kboltz
else
  t=0.d0
end if
write(50) t
write(50) nspinor
write(50) nsymkpt
do isym=1,nsymkpt
  write(50) symkpt(:,:,isym)
end do
write(50) ng_vec
write(50) vgl
write(50) nk
write(50) nstsv
write(50) gkmax
write(50) ngkmax
do ik=1,nkpt
  if (usek(ik)) then
    write(50) vkl(:,ik)
    write(50) ngk(1,ik)
! construct index from G+k-vectors to Yambo's G-vector array
    do igk=1,ngk(1,ik)
      wf_igk(igk)=idx(igkig(igk,1,ik))
    end do
    write(50) wf_igk
    call getevalsv(vkl(:,ik),evalsv(:,ik))
    write(50) evalsv(:,ik)
    call getwfpw(vkl(:,ik),vgkl(:,:,1,ik),wfpw)
    wfpw(ngk(1,ik)+1:,:,:)=0.d0
    write(50) wfpw
  end if
end do
close(50)
write(*,*)
write(*,'("Info(yambo): wrote Yambo input data to YAMBO.OUT")')
write(*,*)
deallocate(usek,useg,idx,wf_igk,vgl,wfpw)
return
end subroutine

