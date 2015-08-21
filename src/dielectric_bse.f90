
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dielectric_bse
use modmain
use modmpi
use modtest
implicit none
! local variables
integer a1,a2,ik1,jk1
integer i1,j1,ist1,jst1
integer iw,i,j,l
real(8) e,eji,t1,t2
complex(8) eta,zv(3),zt1
character(256) fname
! allocatable arrays
real(8), allocatable :: w(:)
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: sigma(:,:,:)
! allocate local arrays
allocate(w(nwdos))
allocate(pmat(3,nstsv,nstsv))
allocate(sigma(3,3,nwdos))
! generate energy grid (starting from zero)
t1=wdos(2)/dble(nwdos)
do iw=1,nwdos
  w(iw)=t1*dble(iw-1)
end do
! i divided by the complex relaxation time
eta=cmplx(0.d0,swidth,8)
sigma(:,:,:)=0.d0
do a2=1,nmbse
  e=evalbse(a2)
  zv(:)=0.d0
! loop over non-reduced k-points
  do ik1=1,nkptnr
! equivalent reduced k-point
    jk1=ikmap(ivk(1,ik1),ivk(2,ik1),ivk(3,ik1))
! read the momentum matrix elements from file
    call getpmat(vkl(:,ik1),pmat)
    do i1=1,nvbse
      ist1=istbse(i1,ik1)
      do j1=1,ncbse
        jst1=jstbse(j1,ik1)
        a1=ijkbse(i1,j1,ik1)
        eji=evalsv(jst1,jk1)-evalsv(ist1,jk1)
        zt1=(e/eji)*hmlbse(a1,a2)
        zv(:)=zv(:)+zt1*pmat(:,ist1,jst1)
      end do
    end do
  end do
  e=evalbse(a2)
  if (abs(e).gt.1.d-8) then
    do i=1,3
      do j=1,3
        zt1=zv(i)*conjg(zv(j))/e
        sigma(i,j,:)=sigma(i,j,:)+zt1/(w(:)-e+eta)+conjg(zt1)/(w(:)+e+eta)
      end do
    end do
  end if
end do
zt1=zi*occmax*wkptnr/omega
sigma(:,:,:)=zt1*sigma(:,:,:)
! write the dielectric function to file
! loop over tensor components
do l=1,noptcomp
  i=optcomp(1,l)
  j=optcomp(2,l)
  t1=0.d0
  if (i.eq.j) t1=1.d0
  write(fname,'("EPSILON_BSE_",2I1,".OUT")') i,j
  open(60,file=trim(fname),action='WRITE',form='FORMATTED')
  do iw=1,nwdos
    t2=t1-fourpi*aimag(sigma(i,j,iw)/(w(iw)+eta))
    write(60,'(2G18.10)') w(iw),t2
  end do
  write(60,'("     ")')
  do iw=1,nwdos
    t2=fourpi*dble(sigma(i,j,iw)/(w(iw)+eta))
    write(60,'(2G18.10)') w(iw),t2
  end do
  close(60)
end do
write(*,*)
write(*,'("Info(dielectric_bse):")')
write(*,'(" dielectric tensor written to EPSILON_BSE_ij.OUT")')
write(*,'(" for components")')
do l=1,noptcomp
  write(*,'("  i = ",I1,", j = ",I1)') optcomp(1:2,l)
end do
! write sigma to test file
call writetest(185,'BSE optical conductivity',nv=nwdos,tol=2.d-2,zva=sigma)
deallocate(w,pmat,sigma)
return
end subroutine

