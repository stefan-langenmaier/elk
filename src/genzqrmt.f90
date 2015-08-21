
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzqrmt(zqrmt)
use modmain
implicit none
! arguments
complex(8), intent(out) :: zqrmt(lmmaxvr,nrcmtmax,natmtot)
! local variables
integer is,ia,ias
integer nrc,irc,l,m,lm
real(8) vecqc(3),qc,tp(2),x
complex(8) zt1
! automatic arrays
real(8) jl(0:lmaxvr)
complex(8) ylm(lmmaxvr),zflm(lmmaxvr)
! allocatable arrays
complex(8), allocatable :: zfmt(:,:)
allocate(zfmt(lmmaxvr,nrcmtmax))
! q-vector in Cartesian coordinates
call r3mv(bvec,vecql,vecqc)
! q-vector length and (theta, phi) coordinates
call sphcrd(vecqc,qc,tp)
! q-vector spherical harmonics
call genylm(lmaxvr,tp,ylm)
do is=1,nspecies
  nrc=nrcmt(is)
  do irc=1,nrc
    x=qc*rcmt(irc,is)
    call sbessel(lmaxvr,x,jl)
    lm=0
    do l=0,lmaxvr
      zt1=fourpi*zil(l)*jl(l)
      do m=-l,l
        lm=lm+1
        zflm(lm)=zt1*conjg(ylm(lm))
      end do
    end do
! convert to spherical coordinates
    call zgemv('N',lmmaxvr,lmmaxvr,zone,zbshtvr,lmmaxvr,zflm,1,zzero, &
     zfmt(:,irc),1)
  end do
! mutiply by phase factor and store in zqrmt
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    x=dot_product(vecqc(:),atposc(:,ia,is))
    zt1=cmplx(cos(x),sin(x),8)
    zqrmt(:,1:nrc,ias)=zt1*zfmt(:,1:nrc)
  end do
end do
deallocate(zfmt)
return
end subroutine

