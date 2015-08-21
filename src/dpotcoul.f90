
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dpotcoul
use modmain
use modphonon
implicit none
! local variables
integer nr,ir,igq0
complex(8) zrho0
! automatic arrays
real(8) vn(nrmtmax)
! allocatable arrays
complex(8), allocatable :: zfmt(:,:),gzfmt(:,:,:)
nr=nrmt(isph)
! solve the complex Poisson's equation in the muffin-tins
call genzvclmt(nrmt,spnrmax,spr,nrmtmax,drhomt,dvclmt)
! compute the monopole potential for effective nuclear charge
call potnucl(ptnucl,nr,spr(:,isph),spzn(isph),vn)
! calculate the gradient of the monopole potential
allocate(zfmt(4,nrmtmax),gzfmt(4,nrmtmax,3))
do ir=1,nr
  zfmt(1,ir)=vn(ir)/y00
  zfmt(2:4,ir)=0.d0
end do
call gradzfmt(1,nr,spr(:,isph),4,nrmtmax,zfmt,gzfmt)
! subtract gradient component corresponding to the phonon polarisation
do ir=1,nr
  dvclmt(2:4,ir,iasph)=dvclmt(2:4,ir,iasph)-gzfmt(2:4,ir,ipph)
end do
deallocate(zfmt,gzfmt)
! solve Poisson's equation in the entire unit cell
if (iqph.eq.iq0) then
  igq0=1
else
  igq0=0
end if
call zpotcoul(nrmt,spnrmax,spr,igq0,gqc,jlgqr,ylmgq,sfacgq,drhoir,nrmtmax, &
 dvclmt,dvclir,zrho0)
return
end subroutine

