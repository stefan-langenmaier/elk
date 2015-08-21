
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getwfpwh(vpl,wfpwh)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
! local variables
integer isym,lspl,ilspl,lspn
integer is,ia,ja,ias,jas
integer ik,ist,ispn,irc,lm,n
integer recl,lmmaxvr_,nrcmtmax_
integer natmtot_,nspinor_,nstsv_
real(8) vkl_(3),si(3,3)
real(8) v(3),det,th,t1
complex(8) su2(2,2),zt1,zt2
! allocatable arrays
complex(8), allocatable :: zflm1(:,:,:)
complex(8), allocatable :: zflm2(:,:)
! find the equivalent k-point number and crystal symmetry element
call findkpt(vpl,isym,ik)
! find the record length
inquire(iolength=recl) vkl_,lmmaxvr_,nrcmtmax_,natmtot_,nspinor_,nstsv_,wfpwh
!$OMP CRITICAL
open(80,file='WFPWH.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
read(80,rec=ik) vkl_,lmmaxvr_,nrcmtmax_,natmtot_,nspinor_,nstsv_,wfpwh
close(80)
!$OMP END CRITICAL
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getwfpwh): differing vectors for k-point ",I8)') ik
  write(*,'(" current   : ",3G18.10)') vkl(:,ik)
  write(*,'(" WFPWH.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (lmmaxvr.ne.lmmaxvr_) then
  write(*,*)
  write(*,'("Error(getwfpwh): differing lmmaxvr for k-point ",I8)') ik
  write(*,'(" current   : ",I8)') lmmaxvr
  write(*,'(" WFPWH.OUT : ",I8)') lmmaxvr_
  write(*,*)
  stop
end if
if (nrcmtmax.ne.nrcmtmax_) then
  write(*,*)
  write(*,'("Error(getwfpwh): differing nrcmtmax for k-point ",I8)') ik
  write(*,'(" current   : ",I8)') nrcmtmax
  write(*,'(" WFPWH.OUT : ",I8)') nrcmtmax_
  write(*,*)
  stop
end if
if (nspinor.ne.nspinor_) then
  write(*,*)
  write(*,'("Error(getwfpwh): differing nspinor for k-point ",I8)') ik
  write(*,'(" current   : ",I8)') nspinor
  write(*,'(" WFPWH.OUT : ",I8)') nspinor_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getwfpwh): differing nstsv for k-point ",I8)') ik
  write(*,'(" current   : ",I8)') nstsv
  write(*,'(" WFPWH.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
! if p = k then return
t1=abs(vpl(1)-vkl(1,ik))+abs(vpl(2)-vkl(2,ik))+abs(vpl(3)-vkl(3,ik))
if (t1.lt.epslat) return
n=nrcmtmax*nspinor*nstsv
allocate(zflm1(lmmaxvr,n,natmtot))
allocate(zflm2(lmmaxvr,n))
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! the inverse of the spatial symmetry rotates k into p
ilspl=isymlat(lspl)
si(:,:)=dble(symlat(:,:,ilspl))
! rotate k-point by inverse symmetry matrix
call r3mtv(si,vkl(:,ik),v)
! make a copy of the wavefunctions with indices rearranged
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    n=0
    do ist=1,nstsv
      do ispn=1,nspinor
        do irc=1,nrcmt(is)
          n=n+1
          zflm1(:,n,ias)=wfpwh(:,irc,ias,ispn,ist)
        end do
      end do
    end do
  end do
end do
do is=1,nspecies
  n=nrcmt(is)*nspinor*nstsv
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! equivalent atom for this symmetry
    ja=ieqatom(ia,is,isym)
    jas=idxas(ja,is)
! phase factor from translation
    t1=-twopi*dot_product(vkl(:,ik),atposl(:,ja,is))
    zt1=cmplx(cos(t1),sin(t1),8)
    t1=twopi*dot_product(v(:),atposl(:,ia,is))
    zt1=zt1*cmplx(cos(t1),sin(t1),8)
! rotate the wavefunctions (active transformation)
    call rotzflm(symlatc(:,:,lspl),0,lmaxvr,n,lmmaxvr,zflm1(:,:,jas),zflm2)
! multiply with phase factor and store
    n=0
    do ist=1,nstsv
      do ispn=1,nspinor
        do irc=1,nrcmt(is)
          n=n+1
          wfpwh(:,irc,ias,ispn,ist)=zt1*zflm2(:,n)
        end do
      end do
    end do
  end do
end do
! apply spin rotation if required
if (spinpol) then
! index to global spin rotation in lattice point group
  lspn=lspnsymc(isym)
! find the SU(2) representation of the spin rotation matrix
  call rotaxang(epslat,symlatc(:,:,lspn),det,v,th)
  call axangsu2(v,th,su2)
! apply SU(2) matrix to spinor wavefunctions (active transformation)
  do ist=1,nstsv
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do irc=1,nrcmt(is)
          do lm=1,lmmaxvr
            zt1=wfpwh(lm,irc,ias,1,ist)
            zt2=wfpwh(lm,irc,ias,2,ist)
            wfpwh(lm,irc,ias,1,ist)=su2(1,1)*zt1+su2(1,2)*zt2
            wfpwh(lm,irc,ias,2,ist)=su2(2,1)*zt1+su2(2,2)*zt2
          end do
        end do
      end do
    end do
  end do
end if
deallocate(zflm1,zflm2)
return
end subroutine

