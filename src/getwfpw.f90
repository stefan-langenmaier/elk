
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getwfpw(vpl,vgpl,wfpw)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vgpl(3,ngkmax,nspnfv)
complex(8), intent(out) :: wfpw(ngkmax,nspinor,nstsv)
! local variables
integer isym,lspl,ilspl,lspn
integer ik,ist,igk,igp,jgp,ig
integer ispn,jspn0,jspn1,i
integer recl,ngkmax_,nspinor_,nstsv_
real(8) vkl_(3),si(3,3)
real(8) v(3),det,th,t1
complex(8) su2(2,2),zt1,zt2
! allocatable arrays
complex(8), allocatable :: wfpwt(:,:,:)
! find the equivalent k-point number and crystal symmetry element
call findkpt(vpl,isym,ik)
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! find the record length
inquire(iolength=recl) vkl_,ngkmax_,nspinor_,nstsv_,wfpw
!$OMP CRITICAL
open(80,file='WFPW.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
read(80,rec=ik) vkl_,ngkmax_,nspinor_,nstsv_,wfpw
close(80)
!$OMP END CRITICAL
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getwfpw): differing vectors for k-point ",I8)') ik
  write(*,'(" current  : ",3G18.10)') vkl(:,ik)
  write(*,'(" WFPW.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (ngkmax.ne.ngkmax_) then
  write(*,*)
  write(*,'("Error(getwfpw): differing nmatmax for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') ngkmax
  write(*,'(" WFPW.OUT : ",I8)') ngkmax_
  write(*,*)
  stop
end if
if (nspinor.ne.nspinor_) then
  write(*,*)
  write(*,'("Error(getwfpw): differing nspinor for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') nspinor
  write(*,'(" WFPW.OUT : ",I8)') nspinor_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getwfpw): differing nstsv for k-point ",I8)') ik
  write(*,'(" current  : ",I8)') nstsv
  write(*,'(" WFPW.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
! if p = k then return
t1=abs(vpl(1)-vkl(1,ik))+abs(vpl(2)-vkl(2,ik))+abs(vpl(3)-vkl(3,ik))
if (t1.lt.epslat) return
!--------------------------------------------------------!
!     translate and rotate wavefunction coefficients     !
!--------------------------------------------------------!
! allocate temporary copy of wavefunction
allocate(wfpwt(ngkmax,nspinor,nstsv))
! the inverse of the spatial symmetry rotates k into p
ilspl=isymlat(lspl)
si(:,:)=dble(symlat(:,:,ilspl))
do ispn=1,nspnfv
  if (spinsprl) then
    jspn0=ispn; jspn1=ispn
  else
    jspn0=1; jspn1=nspinor
  end if
! apply translation operation if required
  if (tvzsymc(isym)) then
! translation vector is zero
    do igk=1,ngk(ispn,ik)
      wfpwt(igk,jspn0:jspn1,:)=wfpw(igk,jspn0:jspn1,:)
    end do
  else
! non-zero translation vector gives a phase factor
    do igk=1,ngk(ispn,ik)
      ig=igkig(igk,ispn,ik)
      t1=-twopi*dot_product(dble(ivg(:,ig)),vtlsymc(:,isym))
      zt1=cmplx(cos(t1),sin(t1),8)
      wfpwt(igk,jspn0:jspn1,:)=zt1*wfpw(igk,jspn0:jspn1,:)
    end do
  end if
! apply spatial rotation operation (passive transformation)
  i=1
  do igk=1,ngk(ispn,ik)
    call r3mtv(si,vgkl(:,igk,ispn,ik),v)
    do igp=i,ngk(ispn,ik)
      t1=abs(v(1)-vgpl(1,igp,ispn)) &
        +abs(v(2)-vgpl(2,igp,ispn)) &
        +abs(v(3)-vgpl(3,igp,ispn))
      if (t1.lt.epslat) then
        wfpw(igp,jspn0:jspn1,:)=wfpwt(igk,jspn0:jspn1,:)
        if (igp.eq.i) i=i+1
        goto 10
      end if
    end do
10 continue
  end do
end do
! apply spin rotation if required
if (spinpol) then
! index to global spin rotation in lattice point group
  lspn=lspnsymc(isym)
! if symmetry element is the identity return
  if (lspn.eq.1) return
! find the SU(2) representation of the spin rotation matrix
  call rotaxang(epslat,symlatc(:,:,lspn),det,v,th)
  call axangsu2(v,th,su2)
! apply SU(2) matrix to spinor wavefunctions (active transformation)
  if (spinsprl) then
! spin-spiral case
    wfpw(:,2,:)=0.d0
    i=1
    do igp=1,ngk(1,ik)
      v(:)=vgpl(:,igp,1)-vqlss(:)
      do jgp=i,ngk(2,ik)
        t1=abs(v(1)-vgpl(1,jgp,2)) &
          +abs(v(2)-vgpl(2,jgp,2)) &
          +abs(v(3)-vgpl(3,jgp,2))
        if (t1.lt.epslat) then
          do ist=1,nstsv
            zt1=wfpw(igp,1,ist)
            zt2=wfpw(jgp,2,ist)
            wfpw(igp,1,ist)=su2(1,1)*zt1+su2(1,2)*zt2
            wfpw(jgp,2,ist)=su2(2,1)*zt1+su2(2,2)*zt2
          end do
          if (jgp.eq.i) i=i+1
          goto 20
        end if
      end do
      wfpw(igp,1,:)=0.d0
20 continue
    end do
  else
! normal spin case
    do ist=1,nstsv
      do igp=1,ngk(1,ik)
        zt1=wfpw(igp,1,ist)
        zt2=wfpw(igp,2,ist)
        wfpw(igp,1,ist)=su2(1,1)*zt1+su2(1,2)*zt2
        wfpw(igp,2,ist)=su2(2,1)*zt1+su2(2,2)*zt2
      end do
    end do
  end if
end if
deallocate(wfpwt)
return
end subroutine

