
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genveffig
! !INTERFACE:
subroutine genveffig
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the Fourier transform of the effective potential in the
!   intersitial region. The potential is first multiplied by the characteristic
!   function which zeros it in the muffin-tins. See routine {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ig,ifg
! allocatable arrays
complex(8), allocatable :: zfft(:)
allocate(zfft(ngrtot))
! Fourier transform effective potential to G-space
zfft(:)=veffir(:)
call zfftifc(3,ngrid,-1,zfft)
! zero Fourier components for |G| > 2*gkmax
do ig=1,ngrtot
  if (gc(ig).gt.2.d0*gkmax) then
    ifg=igfft(ig)
    zfft(ifg)=0.d0
  end if
end do
! Fourier transform back to real space
call zfftifc(3,ngrid,1,zfft)
! multiply by characteristic function
zfft(:)=dble(zfft(:))*cfunir(:)
! Fourier transform back to G-space
call zfftifc(3,ngrid,-1,zfft)
! store in global array
do ig=1,ngvec
  ifg=igfft(ig)
  veffig(ig)=zfft(ifg)
end do
deallocate(zfft)
return
end subroutine
!EOC
