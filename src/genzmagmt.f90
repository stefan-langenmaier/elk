
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genzmagmt(is,wfmt11,wfmt12,wfmt21,wfmt22,ld,zmagmt)
use modmain
implicit none
! arguments
integer, intent(in) :: is
complex(8), intent(in) ::  wfmt11(lmmaxvr,nrcmtmax)
complex(8), intent(in) ::  wfmt12(lmmaxvr,nrcmtmax)
complex(8), intent(in) ::  wfmt21(lmmaxvr,nrcmtmax)
complex(8), intent(in) ::  wfmt22(lmmaxvr,nrcmtmax)
integer, intent(in) :: ld
complex(8), intent(out) :: zmagmt(lmmaxvr,nrcmtmax,ld,ndmag)
! local variables
integer nrc,irc,itp
complex(8) zt1,zt2
nrc=nrcmt(is)
! calculate the z-component of magnetisation: up-up - dn-dn
zmagmt(:,1:nrc,1,ndmag)=conjg(wfmt11(:,1:nrc))*wfmt21(:,1:nrc) &
                       -conjg(wfmt12(:,1:nrc))*wfmt22(:,1:nrc)
! non-collinear case
if (ncmag) then
  do irc=1,nrc
    do itp=1,lmmaxvr
! up-dn spin density
      zt1=conjg(wfmt11(itp,irc))*wfmt22(itp,irc)
! dn-up spin density
      zt2=conjg(wfmt12(itp,irc))*wfmt21(itp,irc)
! calculate the x-component: up-dn + dn-up
      zmagmt(itp,irc,1,1)=zt1+zt2
! calculate the y-component: i*(dn-up - up-dn)
      zt1=zt2-zt1
      zmagmt(itp,irc,1,2)=cmplx(-aimag(zt1),dble(zt1),8)
    end do
  end do
end if
return
end subroutine

