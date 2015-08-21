
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genspecies(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
logical core(maxspst),lorb(maxspst)
integer nz,nmax,lmax,nst
integer ne,nrm,nr,ir,i
integer nlo,ist,jst
integer n(maxspst),l(maxspst),k(maxspst)
! default APW band energy
real(8), parameter :: e0=0.15d0
real(8) mass,zn,t1,t2,t3
real(8) rm,rmin,rmax
real(8) occ(maxspst),eval(maxspst)
character(256) symb,name
! allocatable arrays
real(8), allocatable :: r(:),rho(:),vr(:),rwf(:,:,:)
read(fnum,*,err=20) nz
if (nz.le.0) then
  write(*,*)
  write(*,'("Error(genspecies): atomic number negative : ",I8)') nz
  write(*,*)
  stop
end if
read(fnum,*,err=20) symb,name
read(fnum,*,err=20) mass
! convert from 'atomic mass units' to atomic units
mass=mass*amu
read(fnum,*,err=20) rm
read(fnum,*,err=20) nst
if ((nst.le.0).or.(nst.gt.maxspst)) then
  write(*,*)
  write(*,'("Error(genspecies): nst out of range : ",I8)') nst
  write(*,'(" for species ",A)') trim(name)
  write(*,*)
  stop
end if
ne=0
nmax=1
do ist=1,nst
  read(fnum,*,err=20) n(ist),l(ist),k(ist),i
  ne=ne+i
  occ(ist)=i
  nmax=max(nmax,n(ist))
end do
write(*,'("Info(genspecies): running Z = ",I4,", (",A,")")') nz,trim(name)
if (ne.ne.nz) then
  write(*,*)
  write(*,'("Warning(genspecies): atom not neutral, electron number : ",I4)') ne
end if
! nuclear charge in units of e
zn=-dble(nz)
! set the number of radial mesh points proportional to number of nodes
nrm=100*(nmax+1)
! minimum radial mesh point proportional to 1/sqrt(Z)
rmin=2.d-6/sqrt(dble(nz))
! default effective infinity
rmax=100.d0
do i=1,2
! number of points to effective infinity
  t1=log(rm/rmin)
  t2=log(rmax/rmin)
  t3=dble(nrm)*t2/t1
  nr=int(t3)
  allocate(r(nr),rho(nr),vr(nr),rwf(nr,2,nst))
! generate logarithmic radial mesh
  t2=t1/dble(nrm-1)
  do ir=1,nr
    r(ir)=rmin*exp(dble(ir-1)*t2)
  end do
! solve the Kohn-Sham-Dirac equation for the atom
  call atom(sol,.true.,zn,nst,n,l,k,occ,3,0,nprad,nr,r,eval,rho,vr,rwf)
  do ir=nr,1,-1
    if (rho(ir).gt.1.d-20) then
      rmax=1.5d0*r(ir)
      goto 10
    end if
  end do
10 continue
  deallocate(r,rho,vr,rwf)
end do
! find which states belong to core
do ist=1,nst
  if (eval(ist).lt.ecvcut) then
    core(ist)=.true.
  else
    core(ist)=.false.
  end if
end do
! check that the state for same n and l but different k is also core
do ist=1,nst
  if (core(ist)) then
    do jst=1,nst
      if ((n(ist).eq.n(jst)).and.(l(ist).eq.l(jst))) core(jst)=.true.
    end do
  end if
end do
lmax=1
do ist=1,nst
  if (.not.core(ist)) lmax=max(lmax,l(ist))
end do
! determine the local orbitals
nlo=lmax+1
lorb(:)=.false.
do ist=1,nst
  if (.not.core(ist)) then
    if ((l(ist).eq.0).or.(l(ist).lt.k(ist))) then
      if ((eval(ist).lt.esccut).or.(l(ist).ge.2)) then
        lorb(ist)=.true.
        nlo=nlo+1
      end if
    end if
  end if
end do
! write the species file
open(60,file=trim(symb)//'.in',action='WRITE',form='FORMATTED')
write(60,'(" ''",A,"''",T45,": spsymb")') trim(symb)
write(60,'(" ''",A,"''",T45,": spname")') trim(name)
write(60,'(G14.6,T45,": spzn")') zn
write(60,'(G18.10,T45,": spmass")') mass
write(60,'(G14.6,2F10.4,I6,T45,": sprmin, rmt, sprmax, nrmt")') rmin,rm,rmax,nrm
write(60,'(I4,T45,": spnst")') nst
write(60,'(3I4,G14.6,L1,T45,": spn, spl, spk, spocc, spcore")') n(1),l(1), &
 k(1),occ(1),core(1)
do ist=2,nst
  write(60,'(3I4,G14.6,L1)') n(ist),l(ist),k(ist),occ(ist),core(ist)
end do
write(60,'(I4,T45,": apword")') 1
write(60,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') e0,0,.false.
write(60,'(I4,T45,": nlx")') 0
write(60,'(I4,T45,": nlorb")') nlo
do i=0,lmax
  write(60,'(2I4,T45,": lorbl, lorbord")') i,2
  write(60,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') e0,0,.false.
  write(60,'(F8.4,I4,"  ",L1)') e0,1,.false.
end do
do ist=1,nst
  if (lorb(ist)) then
    write(60,'(2I4,T45,": lorbl, lorbord")') l(ist),3
    write(60,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') e0,0,.false.
    write(60,'(F8.4,I4,"  ",L1)') e0,1,.false.
! subtract 0.25 Ha from eigenvalue because findband searches upwards
    write(60,'(F8.4,I4,"  ",L1)') eval(ist)-0.25d0,0,.true.
  end if
end do
close(60)
return
20 continue
write(*,*)
write(*,'("Error(genspecies): error reading species data")')
write(*,*)
stop
end subroutine

