
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine kerxc(fxcmt,fxcir)
use modmain
use modfxcifc
implicit none
! arguments
real(8), intent(out) :: fxcmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(out) :: fxcir(ngrtot)
! local variables
integer n,is,ia,ias
integer idm,nr,ir,i
real(8) t1
real(8), allocatable :: rho(:),rhoup(:),rhodn(:),mag(:,:)
real(8), allocatable :: fxc(:),fxcuu(:),fxcud(:),fxcdd(:)
! number of independent spin components
n=lmmaxvr*nrmtmax
allocate(rho(n),fxc(n))
if (spinpol) then
  allocate(mag(n,3))
  n=max(n,ngrtot)
  allocate(rhoup(n),rhodn(n))
  allocate(fxcuu(n),fxcud(n),fxcdd(n))
end if
!---------------------------!
!     muffin-tin kernel     !
!---------------------------!
do is=1,nspecies
  nr=nrmt(is)
  n=lmmaxvr*nr
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the density in spherical coordinates
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rhomt(:,:,ias), &
     lmmaxvr,0.d0,rho,lmmaxvr)
    if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
! magnetisation in spherical coordinates
      do idm=1,ndmag
        call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rbshtvr,lmmaxvr, &
         magmt(:,:,ias,idm),lmmaxvr,0.d0,mag(:,idm),lmmaxvr)
      end do
      if (ncmag) then
! non-collinear (use Kubler's trick)
        do i=1,n
! compute rhoup=(rho+|m|)/2 and rhodn=(rho-|m|)/2
          t1=sqrt(mag(i,1)**2+mag(i,2)**2+mag(i,3)**2)
          rhoup(i)=0.5d0*(rho(i)+t1)
          rhodn(i)=0.5d0*(rho(i)-t1)
        end do
      else
! collinear
        do i=1,n
! compute rhoup=(rho+m_z)/2 and rhodn=(rho-m_z)/2
          rhoup(i)=0.5d0*(rho(i)+mag(i,1))
          rhodn(i)=0.5d0*(rho(i)-mag(i,1))
        end do
      end if
! compute fxc
      call fxcifc(fxctype,n=n,rhoup=rhoup,rhodn=rhodn,fxcuu=fxcuu,fxcud=fxcud, &
       fxcdd=fxcdd)
! form the scalar quantity dv/drho
      do i=1,n
        fxc(i)=0.25d0*(fxcuu(i)+2.d0*fxcud(i)+fxcdd(i))
      end do
    else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
      call fxcifc(fxctype,n=n,rho=rho,fxc=fxc)
    end if
! convert fxc to spherical harmonics
    call dgemm('N','N',lmmaxvr,nr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,fxc,lmmaxvr, &
     0.d0,fxcmt(:,:,ias),lmmaxvr)
  end do
end do
!-----------------------------!
!     interstitial kernel     !
!-----------------------------!
if (spinpol) then
!------------------------!
!     spin-polarised     !
!------------------------!
  if (ncmag) then
! non-collinear
    do ir=1,ngrtot
      t1=sqrt(magir(ir,1)**2+magir(ir,2)**2+magir(ir,3)**2)
      rhoup(ir)=0.5d0*(rhoir(ir)+t1)
      rhodn(ir)=0.5d0*(rhoir(ir)-t1)
    end do
  else
! collinear
    do ir=1,ngrtot
      rhoup(ir)=0.5d0*(rhoir(ir)+magir(ir,1))
      rhodn(ir)=0.5d0*(rhoir(ir)-magir(ir,1))
    end do
  end if
  call fxcifc(fxctype,n=ngrtot,rhoup=rhoup,rhodn=rhodn,fxcuu=fxcuu, &
   fxcud=fxcud,fxcdd=fxcdd)
  do ir=1,ngrtot
    fxcir(ir)=0.25d0*(fxcuu(ir)+2.d0*fxcud(ir)+fxcdd(ir))
  end do
else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
  call fxcifc(fxctype,n=ngrtot,rho=rhoir,fxc=fxcir)
end if
deallocate(rho,fxc)
if (spinpol) then
  deallocate(mag,rhoup,rhodn)
  deallocate(fxcuu,fxcud,fxcdd)
end if
return
end subroutine

