
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potcoul
! !INTERFACE:
subroutine potcoul
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the Coulomb potential of the real charge density stored in the
!   global variables {\tt rhomt} and {\tt rhoir} by solving Poisson's equation.
!   These variables are coverted to complex representations and passed to the
!   routine {\tt zpotcoul}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir
real(8) t1
complex(8) zrho0
! automatic arrays
real(8) vn(nrmtmax)
! allocatable arrays
real(8), allocatable :: jlgr(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
allocate(jlgr(0:lnpsd+1,ngvec,nspecies))
allocate(zrhomt(lmmaxvr,nrmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot))
allocate(zvclir(ngrtot))
! convert real muffin-tin charge density to complex spherical harmonic expansion
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmt(is)
    call rtozflm(lmaxvr,rhomt(:,ir,ias),zrhomt(:,ir,ias))
  end do
end do
! store real interstitial charge density in complex array
zrhoir(:)=rhoir(:)
! compute the required spherical Bessel functions
call genjlgpr(lnpsd+1,gc,jlgr)
! solve the complex Poisson's equation in the muffin-tins
call genzvclmt(nrmt,spnrmax,spr,nrmtmax,zrhomt,zvclmt)
! add the nuclear monopole potentials
t1=1.d0/y00
do is=1,nspecies
  call potnucl(ptnucl,nrmt(is),spr(:,is),spzn(is),vn)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      zvclmt(1,ir,ias)=zvclmt(1,ir,ias)+t1*vn(ir)
    end do
  end do
end do
! solve Poisson's equation in the entire unit cell
call zpotcoul(nrmt,spnrmax,spr,1,gc,jlgr,ylmg,sfacg,zrhoir,nrmtmax,zvclmt, &
 zvclir,zrho0)
! convert complex muffin-tin potential to real spherical harmonic expansion
do ias=1,natmtot
  is=idxis(ias)
  do ir=1,nrmt(is)
    call ztorflm(lmaxvr,zvclmt(:,ir,ias),vclmt(:,ir,ias))
  end do
end do
! store complex interstitial potential in real array
vclir(:)=dble(zvclir(:))
! apply constant electric field if required
if (efieldpol) call potefield
deallocate(jlgr,zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine
!EOC

