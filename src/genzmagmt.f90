
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
integer nrc,nrci,irc
integer lmmax,itp
complex(8) z1,z2
nrc=nrcmt(is)
nrci=nrcmtinr(is)
! calculate the z-component of magnetisation: up-up - dn-dn
do irc=1,nrc
  if (irc.le.nrci) then
    lmmax=lmmaxinr
  else
    lmmax=lmmaxvr
  end if
  zmagmt(1:lmmax,irc,1,ndmag)=conjg(wfmt11(1:lmmax,irc))*wfmt21(1:lmmax,irc) &
                             -conjg(wfmt12(1:lmmax,irc))*wfmt22(1:lmmax,irc)
end do
! non-collinear case
if (ncmag) then
  do irc=1,nrc
    if (irc.le.nrci) then
      lmmax=lmmaxinr
    else
      lmmax=lmmaxvr
    end if
    do itp=1,lmmax
! up-dn spin density
      z1=conjg(wfmt11(itp,irc))*wfmt22(itp,irc)
! dn-up spin density
      z2=conjg(wfmt12(itp,irc))*wfmt21(itp,irc)
! calculate the x-component: up-dn + dn-up
      zmagmt(itp,irc,1,1)=z1+z2
! calculate the y-component: i*(dn-up - up-dn)
      z1=z2-z1
      zmagmt(itp,irc,1,2)=cmplx(-aimag(z1),dble(z1),8)
    end do
  end do
end if
return
end subroutine

