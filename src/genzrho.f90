
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genzrho
! !INTERFACE:
subroutine genzrho(tsh,wfmt1,wfmt2,wfir1,wfir2,zrhomt,zrhoir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh    : .true. if the muffin-tin density is to be in spherical harmonics
!            (in,logical)
!   wfmt1  : muffin-tin part of wavefunction 1 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   wfmt2  : muffin-tin part of wavefunction 2 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   wfir1  : interstitial wavefunction 1 (in,complex(ngrtot))
!   wfir2  : interstitial wavefunction 2 (in,complex(ngrtot))
!   zrhomt : muffin-tin charge density in spherical harmonics/coordinates
!            (out,complex(lmmaxvr,nrcmtmax,natmtot))
!   zrhoir : interstitial charge density (out,complex(ngrtot))
! !DESCRIPTION:
!   Calculates the complex overlap charge density from two input wavefunctions:
!   $$ \rho({\bf r})\equiv\Psi_1^{\dag}({\bf r})\cdot\Psi_2({\bf r}). $$
!   Note that the muffin-tin wavefunctions are provided in spherical coordinates
!   and the returned density is either in terms of spherical harmonic
!   coefficients or spherical coordinates when {\tt tsh} is {\tt .true.} or
!   {\tt .false.}, respectively.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
complex(8), intent(in) ::  wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfir1(ngrtot,nspinor)
complex(8), intent(in) ::  wfir2(ngrtot,nspinor)
complex(8), intent(out) :: zrhomt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(out) :: zrhoir(ngrtot)
! local variables
integer is,ia,ias,nrc,n
! allocatable arrays
complex(8), allocatable :: zfmt(:,:)
if (tsh) allocate(zfmt(lmmaxvr,nrcmtmax))
! muffin-tin part
do is=1,nspecies
  nrc=nrcmt(is)
  n=lmmaxvr*nrc
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    if (tsh) then
      if (spinpol) then
        call zvmul2(n,wfmt1(:,:,ias,1),wfmt2(:,:,ias,1),wfmt1(:,:,ias,2), &
         wfmt2(:,:,ias,2),zfmt)
      else
        call zvmul1(n,wfmt1(:,:,ias,1),wfmt2(:,:,ias,1),zfmt)
      end if
! convert to spherical harmonics
      call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr,zfmt, &
       lmmaxvr,zzero,zrhomt(:,:,ias),lmmaxvr)
    else
      if (spinpol) then
        call zvmul2(n,wfmt1(:,:,ias,1),wfmt2(:,:,ias,1),wfmt1(:,:,ias,2), &
         wfmt2(:,:,ias,2),zrhomt(:,:,ias))
      else
        call zvmul1(n,wfmt1(:,:,ias,1),wfmt2(:,:,ias,1),zrhomt(:,:,ias))
      end if
    end if
  end do
end do
if (tsh) deallocate(zfmt)
! interstitial part
if (spinpol) then
! spin-polarised
  call zvmul2(ngrtot,wfir1(:,1),wfir2(:,1),wfir1(:,2),wfir2(:,2),zrhoir)
else
! spin-unpolarised
  call zvmul1(ngrtot,wfir1(:,1),wfir2(:,1),zrhoir)
end if
return
end subroutine
!EOC

subroutine zvmul1(n,a,b,c)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: a(n),b(n)
complex(8), intent(out) :: c(n)
c(:)=conjg(a(:))*b(:)
return
end subroutine

subroutine zvmul2(n,a1,b1,a2,b2,c)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: a1(n),b1(n),a2(n),b2(n)
complex(8), intent(out) :: c(n)
c(:)=conjg(a1(:))*b1(:)+conjg(a2(:))*b2(:)
return
end subroutine

