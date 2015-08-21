
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genephmat(iq,ikp,dvphmt,dvphir,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq,ikp
complex(8), intent(in) :: dvphmt(lmmaxvr,nrcmtmax,natmtot,nbph)
complex(8), intent(in) :: dvphir(ngtot,nbph)
complex(8), intent(out) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer ist,jst,i
real(8) vpql(3)
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:,:),wfirq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
! external functions
complex(8) zfinp
external zfinp
! k+q-vector in lattice coordinates
vpql(:)=vkl(:,ikp)+vql(:,iq)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! generate the wavefunctions for all states at k and k+q
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,nstsv,idx,vkl(:,ikp),wfmt,ngtot,wfir)
allocate(wfmtq(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirq(ngtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,nstsv,idx,vpql,wfmtq,ngtot,wfirq)
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
do ist=1,nstsv
  do jst=1,nstsv
! note that the complex conjugate of the density is found because zfinp
! conjugates the first function
    call genzrho(.true.,.true.,wfmt(:,:,:,:,jst),wfir(:,:,jst), &
     wfmtq(:,:,:,:,ist),wfirq(:,:,ist),zrhomt,zrhoir)
    do i=1,nbph
      ephmat(ist,jst,i)=zfinp(.true.,zrhomt,zrhoir,dvphmt(:,:,:,i),dvphir(:,i))
    end do
  end do
end do
deallocate(wfmt,wfmtq,wfir,wfirq,zrhomt,zrhoir)
return
end subroutine

