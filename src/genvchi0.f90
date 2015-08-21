
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvchi0(ikp,icmp,scsr,vqpl,igq0,gqc,ylmgq,sfacgq,vchi0)
use modmain
implicit none
! local variables
integer, intent(in) :: ikp
integer, intent(in) :: icmp
real(8), intent(in) :: scsr
real(8), intent(in) :: vqpl(3)
integer, intent(in) :: igq0
real(8), intent(in) :: gqc(ngrf)
complex(8), intent(in) :: ylmgq(lmmaxvr,ngrf)
complex(8), intent(in) :: sfacgq(ngrf,natmtot)
complex(8), intent(inout) :: vchi0(nwrf,ngrf,ngrf)
! local variables
integer isym,jkp,jkpq,iw
integer ig,jg,ist,jst
real(8) vkql(3),eij,t1,t2
complex(8) z1,z2
! allocatable arrays
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:,:),wfirq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zw(:),zg(:)
! k+q-vector in lattice coordinates
vkql(:)=vkl(:,ikp)+vqpl(:)
! equivalent reduced k-points for k and k+q
call findkpt(vkl(:,ikp),isym,jkp)
call findkpt(vkql,isym,jkpq)
! generate the wavefunctions for all states at k and k+q
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,.false.,vkl(:,ikp),wfmt,ngtot,wfir)
allocate(wfmtq(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirq(ngtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,.false.,vkql,wfmtq,ngtot,wfirq)
! read the momentum matrix elements from file
allocate(pmat(3,nstsv,nstsv))
call getpmat(vkl(:,ikp),pmat)
! divide by unit cell volume
t1=1.d0/omega
pmat(:,:,:)=t1*pmat(:,:,:)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zw,zg) &
!$OMP PRIVATE(jst,t1,t2,eij) &
!$OMP PRIVATE(iw,ig,jg,z1,z2)
!$OMP DO
do ist=1,nstsv
  if (abs(evalsv(ist,jkp)-efermi).gt.emaxrf) cycle
  allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
  allocate(zw(nwrf),zg(ngrf))
  do jst=1,nstsv
    if (abs(evalsv(jst,jkpq)-efermi).gt.emaxrf) cycle
    t1=wkptnr*omega*(occsv(ist,jkp)-occsv(jst,jkpq))
    if (abs(t1).lt.1.d-8) cycle
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
    do iw=1,nwrf
      zw(iw)=t1/(eij+wrf(iw))
    end do
! compute the complex density in G+q-space
    call genzrho(.true.,.true.,wfmt(:,:,:,:,ist),wfmtq(:,:,:,:,jst), &
     wfir(:,:,ist),wfirq(:,:,jst),zrhomt,zrhoir)
    call zftzf(ngrf,gqc,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zg)
!$OMP CRITICAL
!------------------------!
!     body of matrix     !
!------------------------!
    do jg=1,ngrf
      z1=conjg(zg(jg))
      do ig=1,ngrf
        t1=gqc(ig)*gqc(jg)
        if (t1.gt.1.d-8) then
          z2=(fourpi/t1)*zg(ig)*z1
          call zaxpy(nwrf,z2,zw,1,vchi0(:,ig,jg),1)
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
      vchi0(:,igq0,igq0)=vchi0(:,igq0,igq0)+t1*zw(:)
!-------------------------!
!     wings of matrix     !
!-------------------------!
      t1=-fourpi/eij
      if (icmp.eq.0) then
        z1=(t1/3.d0)*(pmat(1,ist,jst)+pmat(2,ist,jst)+pmat(3,ist,jst))
      else
        z1=t1*pmat(icmp,ist,jst)
      end if
! G = q = 0
      do ig=2,ngrf
        z2=zg(ig)*conjg(z1)/gqc(ig)
        call zaxpy(nwrf,z2,zw,1,vchi0(:,ig,igq0),1)
      end do
! G' = q = 0
      do jg=2,ngrf
        z2=z1*conjg(zg(jg))/gqc(jg)
        call zaxpy(nwrf,z2,zw,1,vchi0(:,igq0,jg),1)
      end do
    end if
!$OMP END CRITICAL
! end loop over jst
  end do
  deallocate(zrhomt,zrhoir,zw,zg)
! end loop over ist
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(pmat,wfmt,wfmtq,wfir,wfirq)
return
end subroutine

