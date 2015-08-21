
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hflocal(vmt,vir,bmt,bir)
use modmain
implicit none
! arguments
real(8), intent(out) :: vmt(lmmaxvr,nrcmtmax,natmtot)
real(8), intent(out) :: vir(ngtot)
real(8), intent(out) :: bmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
real(8), intent(out) :: bir(ngtot,ndmag)
! local variables
integer ld,is,ias
integer idm,ir,irc
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
! compute the Coulomb potential
call potcoul
! generate the exchange-correlation potentials for hybrids
if (hybrid) call potxc
! convert to spherical coordinates and store in output arrays
ld=lmmaxvr*lradstp
if (hybrid) then
! hybrid functional case
  allocate(rfmt(lmmaxvr,nrcmtmax))
  do ias=1,natmtot
    is=idxis(ias)
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      rfmt(:,irc)=vclmt(:,ir,ias)+vxcmt(:,ir,ias)
    end do
    call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr,rfmt, &
     lmmaxvr,0.d0,vmt(:,:,ias),lmmaxvr)
  end do
  deallocate(rfmt)
  vir(:)=(vclir(:)+vxcir(:))*cfunir(:)
  if (spinpol) then
    do idm=1,ndmag
      do ias=1,natmtot
        is=idxis(ias)
        call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
         bxcmt(:,:,ias,idm),ld,0.d0,bmt(:,:,ias,idm),lmmaxvr)
      end do
      bir(:,idm)=bxcir(:,idm)*cfunir(:)
    end do
  end if
else
! normal Hartree-Fock case
  do ias=1,natmtot
    is=idxis(ias)
    call dgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
     vclmt(:,:,ias),ld,0.d0,vmt(:,:,ias),lmmaxvr)
  end do
  vir(:)=vclir(:)*cfunir(:)
end if
return
end subroutine

