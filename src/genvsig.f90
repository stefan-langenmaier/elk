
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvsig
! !INTERFACE:
subroutine genvsig
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the Fourier transform of the Kohn-Sham effective potential in the
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
real(8) gm2
! allocatable arrays
complex(8), allocatable :: zfft(:)
allocate(zfft(ngtot))
if (trimvg) then
! trim the Fourier components of vsir for |G| > gmaxvr/2
  zfft(:)=vsir(:)
  call zfftifc(3,ngridg,-1,zfft)
  gm2=gmaxvr/2.d0
  do ig=1,ngtot
    if (gc(ig).gt.gm2) then
      ifg=igfft(ig)
      zfft(ifg)=0.d0
    end if
  end do
! Fourier transform back to real-space
  call zfftifc(3,ngridg,1,zfft)
! multiply trimmed potential by characteristic function in real-space
  zfft(:)=dble(zfft(:))*cfunir(:)
else
! multiply potential by characteristic function in real-space
  zfft(:)=vsir(:)*cfunir(:)
end if
! Fourier transform to G-space
call zfftifc(3,ngridg,-1,zfft)
! store in global array
do ig=1,ngvec
  ifg=igfft(ig)
  vsig(ig)=zfft(ifg)
end do
deallocate(zfft)
return
end subroutine
!EOC

