
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genepsqk(iq,ikp,gqc,expqmt,eps)
use modmain
implicit none
! local variables
integer, intent(in) :: iq
integer, intent(in) :: ikp
real(8), intent(in) :: gqc(ngrpa)
complex(8), intent(in) :: expqmt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(inout) :: eps(nwrpa,ngrpa,ngrpa)
! local variables
integer ispn,is,ia,ias,irc
integer ist,jst,iw,ig,jg
integer isym,jkp,jkpq
real(8) vpql(3),eij,t1
complex(8) zt1,zt2
! allocatable arrays
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: wfmtp(:,:,:,:,:)
complex(8), allocatable :: wfmtpq(:,:,:,:,:)
complex(8), allocatable :: wfirp(:,:,:)
complex(8), allocatable :: wfirpq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zw(:)
complex(8), allocatable :: zv(:)
! external functions
complex(8) zfinp
external zfinp
! p+q vector in lattice coordinates
vpql(:)=vkl(:,ikp)+vql(:,iq)
! equivalent reduced k-points for p and p+q
call findkpt(vkl(:,ikp),isym,jkp)
call findkpt(vpql,isym,jkpq)
! generate the wavefunctions for all states at p and p+q
allocate(wfmtp(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirp(ngrtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,vkl(:,ikp),wfmtp,ngrtot,wfirp)
allocate(wfmtpq(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirpq(ngrtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,vpql,wfmtpq,ngrtot,wfirpq)
! read the momentum matrix elements from file
allocate(pmat(3,nstsv,nstsv))
call getpmat(vkl(:,ikp),pmat)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zw,zv) &
!$OMP PRIVATE(ispn,is,ia,ias,irc) &
!$OMP PRIVATE(jst,t1,eij,iw,ig,jg) &
!$OMP PRIVATE(zt1,zt2)
!$OMP DO
do ist=1,nstsv
  allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
  allocate(zrhoir(ngrtot))
  allocate(zw(nwrpa))
  allocate(zv(ngrpa))
! multiply wfmtp with exp(iq.r)
  do ispn=1,nspinor
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do irc=1,nrcmt(is)
          wfmtp(:,irc,ias,ispn,ist)=wfmtp(:,irc,ias,ispn,ist)*expqmt(:,irc,ias)
        end do
      end do
    end do
  end do
  do jst=1,nstsv
    t1=(wkptnr/omega)*occsv(ist,jkp)*(1.d0-occsv(jst,jkpq)/occmax)
    if (t1.gt.1.d-8) then
      eij=evalsv(ist,jkp)-(evalsv(jst,jkpq)+scissor)
! frequency dependent part in response function formula for all RPA frequencies
! (note: this formula is unsuitable for systems without time-reversal symmetry)
      do iw=1,nwrpa
        zw(iw)=t1*(1.d0/(eij+wrpa(iw))+1.d0/(eij-wrpa(iw)))
      end do
! compute the complex density, note complex conjugate is computed because zfinp
! conjugates the first function
      call genzrho(.false.,wfmtp(:,:,:,:,ist),wfmtpq(:,:,:,:,jst), &
       wfirp(:,:,ist),wfirpq(:,:,jst),zrhomt,zrhoir)
      do ig=1,ngrpa
        zv(ig)=zfinp(.false.,zrhomt,expgmt(:,:,:,ig),zrhoir,expgir(:,ig))
      end do
! q non-zero case
      if (iq.ne.iq0) then
! add to the dielectric function (in general epsilon is not Hermitian)
!$OMP CRITICAL
        do ig=1,ngrpa
          zt1=conjg(zv(ig))
          do jg=1,ngrpa
            t1=-fourpi/(gqc(ig)*gqc(jg))
            zt2=t1*zt1*zv(jg)
            eps(:,ig,jg)=eps(:,ig,jg)+zt2*zw(:)
          end do
        end do
!$OMP END CRITICAL
      end if
! special case of q = 0
      if (iq.eq.iq0) then
!$OMP CRITICAL
! head of matrix: G = G' = q = 0
        t1=sum(dble(pmat(:,ist,jst))**2+aimag(pmat(:,ist,jst))**2)
        t1=-fourpi*t1/(3.d0*eij**2)
        eps(:,1,1)=eps(:,1,1)+t1*zw(:)
! wings of matrix
        t1=fourpi/(3.d0*eij)
        zt1=t1*(pmat(1,ist,jst)+pmat(2,ist,jst)+pmat(3,ist,jst))
! G = q = 0
        do ig=2,ngrpa
          zt2=zt1*conjg(zv(ig))/gqc(ig)
          eps(:,ig,1)=eps(:,ig,1)+zt2*zw(:)
        end do
! G' = q = 0
        do jg=2,ngrpa
          zt2=conjg(zt1)*zv(jg)/gqc(jg)
          eps(:,1,jg)=eps(:,1,jg)+zt2*zw(:)
        end do
!$OMP END CRITICAL
      end if
    end if
! end loop over jst
  end do
  deallocate(zrhomt,zrhoir,zw,zv)
! end loop over ist
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(pmat,wfmtp,wfmtpq,wfirp,wfirpq)
return
end subroutine

