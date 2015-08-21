
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genexpmt(vpc,expmt)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpc(3)
complex(8), intent(out) :: expmt(lmmaxvr,nrcmtmax,natmtot)
! local variables
integer is,ia,ias
integer nrc,irc,l,m,lm
real(8) pc,tp(2),x
complex(8) zt1
! automatic arrays
real(8) jl(0:lmaxvr)
complex(8) ylm(lmmaxvr)
! allocatable arrays
complex(8), allocatable :: zfmt1(:,:),zfmt2(:,:)
allocate(zfmt1(lmmaxvr,nrcmtmax),zfmt2(lmmaxvr,nrcmtmax))
! p-vector length and (theta, phi) coordinates
call sphcrd(vpc,pc,tp)
! p-vector spherical harmonics
call genylm(lmaxvr,tp,ylm)
do is=1,nspecies
  nrc=nrcmt(is)
  do irc=1,nrc
    x=pc*rcmt(irc,is)
    call sbessel(lmaxvr,x,jl)
    lm=0
    do l=0,lmaxvr
      zt1=fourpi*jl(l)*zil(l)
      do m=-l,l
        lm=lm+1
        zfmt1(lm,irc)=zt1*conjg(ylm(lm))
      end do
    end do
  end do
! convert to spherical harmonics
  call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zbshtvr,lmmaxvr,zfmt1,lmmaxvr, &
   zzero,zfmt2,lmmaxvr)
! mutiply by phase factors and store for all atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    x=dot_product(vpc(:),atposc(:,ia,is))
    zt1=cmplx(cos(x),sin(x),8)
    expmt(:,1:nrc,ias)=zt1*zfmt2(:,1:nrc)
  end do
end do
deallocate(zfmt1,zfmt2)
return
end subroutine

