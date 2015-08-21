
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: seceqn
subroutine seceqn(ik,evalfv,evecfv,evecsv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (out,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Solves the first- and second-variational secular equations. See routines
!   {\tt match}, {\tt seceqnfv}, {\tt seceqnss} and {\tt seceqnsv}.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(out) :: evalfv(nstfv,nspnfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer jspn
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
! loop over first-variational spins (nspnfv=2 for spin-spirals only)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do jspn=1,nspnfv
! find the matching coefficients
  call match(ngk(jspn,ik),gkc(:,jspn,ik),tpgkc(:,:,jspn,ik), &
   sfacgk(:,:,jspn,ik),apwalm(:,:,:,:,jspn))
! solve the first-variational secular equation
  if (tseqit) then
! iteratively
    call seceqnit(nmat(jspn,ik),ngk(jspn,ik),igkig(:,jspn,ik),vkl(:,ik), &
     vgkl(:,:,jspn,ik),vgkc(:,:,jspn,ik),apwalm(:,:,:,:,jspn),evalfv(:,jspn), &
     evecfv(:,:,jspn))
  else
! directly
    call seceqnfv(nmat(jspn,ik),ngk(jspn,ik),igkig(:,jspn,ik),vkc(:,ik), &
     vgkc(:,:,jspn,ik),apwalm(:,:,:,:,jspn),evalfv(:,jspn),evecfv(:,:,jspn))
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
if (spinsprl) then
! solve the spin-spiral second-variational secular equation
  call seceqnss(ngk(:,ik),igkig(:,:,ik),apwalm,evalfv,evecfv,evalsv(:,ik), &
   evecsv)
else
! solve the second-variational secular equation
  call seceqnsv(ngk(1,ik),igkig(:,1,ik),vgkc(:,:,1,ik),apwalm,evalfv,evecfv, &
   evalsv(:,ik),evecsv)
end if
deallocate(apwalm)
return
end subroutine
!EOC

