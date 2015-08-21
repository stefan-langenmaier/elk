
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvchi0(ikp,icmp,scsr,vqpl,gqc,expqmt,vchi0)
use modmain
implicit none
! local variables
integer, intent(in) :: ikp
integer, intent(in) :: icmp
real(8), intent(in) :: scsr
real(8), intent(in) :: vqpl(3)
real(8), intent(in) :: gqc(ngrpa)
complex(8), intent(in) :: expqmt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(inout) :: vchi0(ngrpa,ngrpa,nwrpa)
! local variables
integer ispn,ias,is,irc
integer ist,jst,iw,ig,jg
integer isym,jkp,jkpq,igq0
real(8) vkql(3),eij,t1,t2
complex(8) zt1,zt2
! allocatable arrays
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: wfmtp(:,:,:,:,:),wfmtpq(:,:,:,:,:)
complex(8), allocatable :: wfirp(:,:,:),wfirpq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zw(:),zg(:)
! external functions
complex(8) zfinp
external zfinp
! k+q-vector in lattice coordinates
vkql(:)=vkl(:,ikp)+vqpl(:)
! equivalent reduced k-points for k and k+q
call findkpt(vkl(:,ikp),isym,jkp)
call findkpt(vkql,isym,jkpq)
! generate the wavefunctions for all states at k and k+q
allocate(wfmtp(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirp(ngrtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,vkl(:,ikp),wfmtp,ngrtot,wfirp)
allocate(wfmtpq(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirpq(ngrtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,vkql,wfmtpq,ngrtot,wfirpq)
! read the momentum matrix elements from file
allocate(pmat(3,nstsv,nstsv))
call getpmat(vkl(:,ikp),pmat)
! find the shortest G+q-vector
call findigp0(ngrpa,gqc,igq0)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zw,zg) &
!$OMP PRIVATE(ispn,ias,is,irc,jst) &
!$OMP PRIVATE(t1,t2,eij,iw,ig,jg) &
!$OMP PRIVATE(zt1,zt2)
!$OMP DO
do ist=1,nstsv
  allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngrtot))
  allocate(zw(nwrpa),zg(ngrpa))
! multiply wfmtp with exp(iq.r)
  do ispn=1,nspinor
    do ias=1,natmtot
      is=idxis(ias)
      do irc=1,nrcmt(is)
        wfmtp(:,irc,ias,ispn,ist)=wfmtp(:,irc,ias,ispn,ist)*expqmt(:,irc,ias)
      end do
    end do
  end do
  do jst=1,nstsv
    t1=(wkptnr/omega)*(occsv(ist,jkp)-occsv(jst,jkpq))
    if (abs(t1).gt.1.d-8) then
      eij=evalsv(ist,jkp)-evalsv(jst,jkpq)
! scissor operator
      if (abs(scsr).gt.1.d-8) then
        t2=eij
        if (eij.gt.0.d0) then
          eij=eij+scsr
        else
          eij=eij-scsr
        end if
        t2=eij/t2
! scale the momentum matrix elements
        pmat(:,ist,jst)=t2*pmat(:,ist,jst)
      end if
! frequency-dependent part in response function formula for all frequencies
      do iw=1,nwrpa
        zw(iw)=t1/(eij+wrpa(iw))
      end do
! compute the complex density in G+q-space
      call genzrho(.false.,spinpol,wfmtp(:,:,:,:,ist),wfmtpq(:,:,:,:,jst), &
       wfirp(:,:,ist),wfirpq(:,:,jst),zrhomt,zrhoir)
      do ig=1,ngrpa
        zg(ig)=zfinp(.false.,zrhomt,expgmt(:,:,:,ig),zrhoir,expgir(:,ig))
      end do
!$OMP CRITICAL
!------------------------!
!     body of matrix     !
!------------------------!
      do ig=1,ngrpa
        zt1=conjg(zg(ig))
        do jg=1,ngrpa
          t1=gqc(ig)*gqc(jg)
          if (t1.gt.1.d-8) then
            zt2=(fourpi/t1)*zt1*zg(jg)
            vchi0(ig,jg,:)=vchi0(ig,jg,:)+zt2*zw(:)
          end if
        end do
      end do
! special case of q = 0
      if ((gqc(igq0).lt.epslat).and.(abs(eij).gt.1.d-8)) then
!----------------------------------------!
!     head of matrix: G = G' = q = 0     !
!----------------------------------------!
        if (icmp.eq.0) then
! trace of dielectric tensor
          t1=sum(dble(pmat(:,ist,jst))**2+aimag(pmat(:,ist,jst))**2)/3.d0
        else
! particular macroscopic component
          t1=dble(pmat(icmp,ist,jst))**2+aimag(pmat(icmp,ist,jst))**2
        end if
        t1=fourpi*t1/eij**2
        vchi0(igq0,igq0,:)=vchi0(igq0,igq0,:)+t1*zw(:)
!-------------------------!
!     wings of matrix     !
!-------------------------!
        t1=-fourpi/eij
        if (icmp.eq.0) then
          zt1=(t1/3.d0)*(pmat(1,ist,jst)+pmat(2,ist,jst)+pmat(3,ist,jst))
        else
          zt1=t1*pmat(icmp,ist,jst)
        end if
! G = q = 0
        do ig=2,ngrpa
          zt2=conjg(zg(ig))*zt1/gqc(ig)
          vchi0(ig,igq0,:)=vchi0(ig,igq0,:)+zt2*zw(:)
        end do
! G' = q = 0
        do jg=2,ngrpa
          zt2=conjg(zt1)*zg(jg)/gqc(jg)
          vchi0(igq0,jg,:)=vchi0(igq0,jg,:)+zt2*zw(:)
        end do
      end if
!$OMP END CRITICAL
    end if
! end loop over jst
  end do
  deallocate(zrhomt,zrhoir,zw,zg)
! end loop over ist
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(pmat,wfmtp,wfmtpq,wfirp,wfirpq)
return
end subroutine

