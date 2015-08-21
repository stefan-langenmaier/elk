
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genephmat(iq,ik,dvphmt,dvphir,ephmat)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: iq,ik
complex(8), intent(in) :: dvphmt(lmmaxvr,nrcmtmax,natmtot,nbph)
complex(8), intent(in) :: dvphir(ngtot,nbph)
complex(8), intent(out) :: ephmat(nstsv,nstsv,nbph)
! local variables
integer jk,jkq,isym
integer nst,nstq,i
integer ist,jst,kst,lst
real(8), parameter :: eps=1.d-6
real(8) vpql(3),t1
! automatic arrays
integer idx(nstsv),idxq(nstsv)
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:,:),wfirq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
! external functions
complex(8) zfinp
external zfinp
! equivalent reduced k-point
jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! k+q-vector in lattice coordinates
vpql(:)=vkl(:,ik)+vql(:,iq)
! find reduced k-point index corresponding to k+q
call findkpt(vpql,isym,jkq)
! index to states near Fermi energy
nst=0
nstq=0
do ist=1,nstsv
  t1=abs(occsv(ist,jk))/occmax
  if ((t1.gt.eps).and.(t1.lt.1.d0-eps)) then
    nst=nst+1
    idx(nst)=ist
  end if
  t1=abs(occsv(ist,jkq))/occmax
  if ((t1.gt.eps).and.(t1.lt.1.d0-eps)) then
    nstq=nstq+1
    idxq(nstq)=ist
  end if
end do
! generate the wavefunctions for all states at k and k+q
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nst))
allocate(wfir(ngtot,nspinor,nst))
call genwfsvp(.false.,.false.,nst,idx,vkl(:,ik),wfmt,ngtot,wfir)
allocate(wfmtq(lmmaxvr,nrcmtmax,natmtot,nspinor,nstq))
allocate(wfirq(ngtot,nspinor,nstq))
call genwfsvp(.false.,.false.,nstq,idxq,vpql,wfmtq,ngtot,wfirq)
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
ephmat(:,:,:)=0.d0
do ist=1,nstq
  kst=idxq(ist)
  do jst=1,nst
    lst=idx(jst)
! note that the complex conjugate of the density is found because zfinp
! conjugates the first function
    call genzrho(.true.,.true.,wfmt(:,:,:,:,jst),wfir(:,:,jst), &
     wfmtq(:,:,:,:,ist),wfirq(:,:,ist),zrhomt,zrhoir)
    do i=1,nbph
      ephmat(kst,lst,i)=zfinp(zrhomt,zrhoir,dvphmt(:,:,:,i),dvphir(:,i))
    end do
  end do
end do
deallocate(wfmt,wfmtq,wfir,wfirq,zrhomt,zrhoir)
return
end subroutine

