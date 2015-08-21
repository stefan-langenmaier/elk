
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genephmat(iq,ikp,dvphmt,dvphir,ephmat)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
integer, intent(in) :: ikp
complex(8), intent(in) :: dvphmt(lmmaxvr,nrcmtmax,natmtot,3*natmtot)
complex(8), intent(in) :: dvphir(ngtot,3*natmtot)
complex(8), intent(out) :: ephmat(nstsv,nstsv,3*natmtot)
! local variables
integer nb,ist,jst,i
real(8) vkql(3)
! allocatable arrays
complex(8), allocatable :: wfmtp(:,:,:,:,:),wfmtpq(:,:,:,:,:)
complex(8), allocatable :: wfirp(:,:,:),wfirpq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
! external functions
complex(8) zfinp
external zfinp
! number of phonon branches
nb=3*natmtot
! k+q-vector in lattice coordinates
vkql(:)=vkl(:,ikp)+vql(:,iq)
! generate the wavefunctions for all states at k and k+q
allocate(wfmtp(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirp(ngtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,.false.,vkl(:,ikp),wfmtp,ngtot,wfirp)
allocate(wfmtpq(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirpq(ngtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,.false.,vkql,wfmtpq,ngtot,wfirpq)
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
do ist=1,nstsv
  do jst=1,nstsv
! note that the complex conjugate of the density is found because zfinp
! conjugates the first function
    call genzrho(.true.,spinpol,wfmtp(:,:,:,:,jst),wfmtpq(:,:,:,:,ist), &
     wfirp(:,:,jst),wfirpq(:,:,ist),zrhomt,zrhoir)
    do i=1,nb
      ephmat(ist,jst,i)=zfinp(.true.,zrhomt,dvphmt(:,:,:,i),zrhoir,dvphir(:,i))
    end do
  end do
end do
deallocate(wfmtp,wfmtpq,wfirp,wfirpq,zrhomt,zrhoir)
return
end subroutine

