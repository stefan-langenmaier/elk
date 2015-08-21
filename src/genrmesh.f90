
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genrmesh
! !INTERFACE:
subroutine genrmesh
! !USES:
use modmain
use modvars
! !DESCRIPTION:
!   Generates the coarse and fine radial meshes for each atomic species in the
!   crystal. Also determines which points are in the inner part of the
!   muffin-tin using the value of {\tt fracinr}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ir,irc
real(8) t1,t2
! estimate the number of radial mesh points to infinity
spnrmax=1
do is=1,nspecies
! logarithmic mesh
  t1=log(sprmax(is)/sprmin(is))/log(rmt(is)/sprmin(is))
  t2=dble(nrmt(is)-1)*t1
  spnr(is)=nint(t2)+1
  spnrmax=max(spnrmax,spnr(is))
end do
! allocate the global radial mesh arrays
if (allocated(spr)) deallocate(spr)
allocate(spr(spnrmax,nspecies))
if (allocated(rcmt)) deallocate(rcmt)
allocate(rcmt(nrcmtmax,nspecies))
! generate the radial meshes
do is=1,nspecies
  t1=1.d0/dble(nrmt(is)-1)
! logarithmic mesh
  t2=log(rmt(is)/sprmin(is))
  do ir=1,spnr(is)
    spr(ir,is)=sprmin(is)*exp(dble(ir-1)*t1*t2)
  end do
end do
! set up the coarse radial meshes and find the inner part of the muffin-tin
! where rho is calculated with lmaxinr
do is=1,nspecies
  t1=fracinr*rmt(is)
  nrmtinr(is)=0
  nrcmtinr(is)=0
  irc=0
  do ir=1,nrmt(is),lradstp
    irc=irc+1
    rcmt(irc,is)=spr(ir,is)
    if (spr(ir,is).lt.t1) then
      nrmtinr(is)=ir
      nrcmtinr(is)=irc
    end if
  end do
end do
! determine at which radius the Hellmann-Feynman force is to be evaluated
do is=1,nspecies
  do ir=1,nrmt(is)
    if (spr(ir,is).lt.radfhf) irfhf(is)=ir
  end do
end do
! write to VARIABLES.OUT
call writevars('spnr',nv=nspecies,iva=spnr)
call writevars('nrmt',nv=nspecies,iva=nrmt)
call writevars('nrmtinr',nv=nspecies,iva=nrmtinr)
call writevars('lradstp',iv=lradstp)
call writevars('nrcmt',nv=nspecies,iva=nrcmt)
call writevars('nrcmtinr',nv=nspecies,iva=nrcmtinr)
do is=1,nspecies
  call writevars('spr',nv=nrmt(is),rva=spr(:,is))
end do
return
end subroutine
!EOC
