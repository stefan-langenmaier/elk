
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine energynn
use modmain
implicit none
! local variables
integer is,ia,ias,ir
real(8) t1
complex(8) zrho0
! automatic arrays
real(8) vn(nrmtmax),vn0(nspecies)
! allocatable arrays
real(8), allocatable :: jlgr(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
allocate(jlgr(0:lnpsd+1,ngvec,nspecies))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot))
allocate(zvclir(ngrtot))
! set the density to zero
zrhoir(:)=0.d0
! compute the required spherical Bessel functions
call genjlgpr(lnpsd+1,gc,jlgr)
! generate the nuclear monopole potentials
t1=1.d0/y00
do is=1,nspecies
  call potnucl(ptnucl,nrmt(is),spr(:,is),spzn(is),vn)
  vn0(is)=vn(1)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      zvclmt(1,ir,ias)=t1*vn(ir)
      zvclmt(2:,ir,ias)=0.d0
    end do
  end do
end do
! solve the complex Poisson's equation
call zpotcoul(nrmt,spnrmax,spr,1,gc,jlgr,ylmg,sfacg,zrhoir,nrmtmax,zvclmt, &
 zvclir,zrho0)
! compute the nuclear-nuclear energy
engynn=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    t1=dble(zvclmt(1,1,ias))*y00-vn0(is)
    engynn=engynn+spzn(is)*t1
  end do
end do
engynn=0.5d0*engynn
deallocate(jlgr,zrhoir,zvclmt,zvclir)
return
end subroutine

