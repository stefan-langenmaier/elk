
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genphsc(p,dph)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: p
real(8), intent(in) :: dph
! local variables
integer js,ja,na,i,n,iv(3)
integer i1,i2,i3,scl(3,3)
real(8) v1(3),v2(3),dmin,t1
if ((p.ne.0).and.(p.ne.1)) then
  write(*,*)
  write(*,'("Error(genphsc): phase (p) should be 0 or 1 : ",I8)') p
  write(*,*)
  stop
end if
! check for Gamma-point phonon
if ((ivq(1,iqph).eq.0).and.(ivq(2,iqph).eq.0).and.(ivq(3,iqph).eq.0)) then
  scl(:,:)=0
  scl(1,1)=1
  scl(2,2)=1
  scl(3,3)=1
  nphsc=1
  goto 10
end if
! find the first lattice vector
dmin=1.d8
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      t1=dble(i1)*vql(1,iqph)+dble(i2)*vql(2,iqph)+dble(i3)*vql(3,iqph)
      if (abs(t1-nint(t1)).lt.epslat) then
        v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
        t1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
        if ((t1.lt.dmin).and.(t1.gt.epslat)) then
          scl(1,1)=i1
          scl(2,1)=i2
          scl(3,1)=i3
          dmin=t1
        end if
      end if
    end do
  end do
end do
! find the second lattice vector
dmin=1.d8
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      t1=dble(i1)*vql(1,iqph)+dble(i2)*vql(2,iqph)+dble(i3)*vql(3,iqph)
      if (abs(t1-nint(t1)).lt.epslat) then
! area defined by first two lattice vectors
        n=(i2*scl(3,1)-i3*scl(2,1))**2 &
         +(i3*scl(1,1)-i1*scl(3,1))**2 &
         +(i1*scl(2,1)-i2*scl(1,1))**2
        if (n.ne.0) then
          v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
          t1=v1(1)**2+v1(2)**2+v1(3)**2
          if (t1.lt.dmin) then
            scl(1,2)=i1
            scl(2,2)=i2
            scl(3,2)=i3
            dmin=t1
          end if
        end if
      end if
    end do
  end do
end do
! find the third lattice vector
nphsc=0
dmin=1.d8
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      t1=dble(i1)*vql(1,iqph)+dble(i2)*vql(2,iqph)+dble(i3)*vql(3,iqph)
      if (abs(t1-nint(t1)).lt.epslat) then
! number of primitive unit cells in supercell
        n=scl(1,2)*(i2*scl(3,1)-i3*scl(2,1)) &
         +scl(2,2)*(i3*scl(1,1)-i1*scl(3,1)) &
         +scl(3,2)*(i1*scl(2,1)-i2*scl(1,1))
        if (n.ne.0) then
          v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
          t1=v1(1)**2+v1(2)**2+v1(3)**2
          if (t1.lt.dmin) then
            nphsc=abs(n)
            scl(1,3)=i1
            scl(2,3)=i2
            scl(3,3)=i3
            dmin=t1
          end if
        end if
      end if
    end do
  end do
end do
if (nphsc.eq.0) then
  write(*,*)
  write(*,'("Error(genphsc): unable to generate supercell")')
  write(*,*)
  stop
end if
10 continue
! new lattice vectors
do i=1,3
  avec(:,i)=dble(scl(1,i))*avec0(:,1) &
           +dble(scl(2,i))*avec0(:,2) &
           +dble(scl(3,i))*avec0(:,3)
end do
! inverse of lattice vector matrix
call r3minv(avec,ainv)
! generate offset vectors for each primitive cell in the supercell
n=1
vphsc(:,1)=0.d0
do i1=-ngridq(1),ngridq(1)
  do i2=-ngridq(2),ngridq(2)
    do i3=-ngridq(3),ngridq(3)
      if (n.eq.nphsc) goto 30
      v1(:)=dble(i1)*avec0(:,1)+dble(i2)*avec0(:,2)+dble(i3)*avec0(:,3)
      call r3mv(ainv,v1,v2)
      call r3frac(epslat,v2,iv)
      call r3mv(avec,v2,v1)
      do i=1,n
        t1=abs(v1(1)-vphsc(1,i))+abs(v1(2)-vphsc(2,i))+abs(v1(3)-vphsc(3,i))
        if (t1.lt.epslat) goto 20
      end do
      n=n+1
      vphsc(:,n)=v1(:)
20 continue
    end do
  end do
end do
write(*,*)
write(*,'("Error(genphsc): unable to generate supercell")')
write(*,*)
stop
30 continue
! set up new atomic positions
do js=1,nspecies
  na=0
  do ja=1,natoms(js)
    do i=1,nphsc
      na=na+1
      if (na.gt.maxatoms) then
        write(*,*)
        write(*,'("Error(genphsc): too many atoms in supercell : ",I8)') na
        write(*,'(" for species ",I4)') js
        write(*,'("Adjust maxatoms in modmain and recompile code")')
        write(*,*)
        stop
      end if
      v1(:)=vphsc(:,i)+atposc0(:,ja,js)
! add small periodic displacement
      if ((isph.eq.js).and.(iaph.eq.ja)) then
        t1=dot_product(vqc(:,iqph),vphsc(:,i))
        if (p.eq.0) then
          v1(ipph)=v1(ipph)+dph*cos(t1)
        else
          v1(ipph)=v1(ipph)+dph*sin(t1)
        end if
      end if
! convert to new lattice coordinates
      call r3mv(ainv,v1,atposl(:,na,js))
      call r3frac(epslat,atposl(:,na,js),iv)
    end do
  end do
  natoms(js)=na
end do
! muffin-tin magnetic fields should be zero
bfcmt(:,:,:)=0.d0
return
end subroutine

