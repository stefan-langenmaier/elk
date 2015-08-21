
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genchi0(iq,ikp,chi0)
use modmain
implicit none
! local variables
integer, intent(in) :: iq
integer, intent(in) :: ikp
complex(8), intent(inout) :: chi0(nwrpa,ngrpa,ngrpa)
! local variables
integer isym,jkp,jkpq
integer ispn,jspn,ist,jst,iw
integer ngp(nspnfv),ngpq(nspnfv)
integer igp,igpq,ifg,ig,jg
real(8) vpql(3),vpqc(3)
real(8) vl(3),vc(3),eij,t0,t1
complex(8) zt1,zt2
! allocatable arrays
integer, allocatable :: igpig(:,:)
integer, allocatable :: igpqig(:,:)
real(8), allocatable :: vgpl(:,:,:)
real(8), allocatable :: vgpql(:,:,:)
real(8), allocatable :: vgpc(:,:,:)
real(8), allocatable :: vgpqc(:,:,:)
complex(8), allocatable :: wfpwp(:,:,:)
complex(8), allocatable :: wfpwpq(:,:,:)
complex(8), allocatable :: wfirp(:,:,:)
complex(8), allocatable :: wfirpq(:,:,:)
complex(8), allocatable :: zw(:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zv(:)
allocate(igpig(ngkmax,nspnfv))
allocate(igpqig(ngkmax,nspnfv))
allocate(vgpl(3,ngkmax,nspnfv))
allocate(vgpql(3,ngkmax,nspnfv))
allocate(vgpc(3,ngkmax,nspnfv))
allocate(vgpqc(3,ngkmax,nspnfv))
allocate(wfpwp(ngkmax,nspinor,nstsv))
allocate(wfpwpq(ngkmax,nspinor,nstsv))
allocate(wfirp(ngrtot,nspinor,nstsv))
allocate(wfirpq(ngrtot,nspinor,nstsv))
! p+q vector in lattice and Cartesian coordinates
vpql(:)=vkl(:,ikp)+vql(:,iq)
vpqc(:)=vkc(:,ikp)+vqc(:,iq)
! equivalent reduced k-points for p and p+q
call findkpt(vkl(:,ikp),isym,jkp)
call findkpt(vpql,isym,jkpq)
! generate the G+p and G+p+q vectors
do ispn=1,nspnfv
  vl(:)=vkl(:,ikp)
  vc(:)=vkc(:,ikp)
  if (spinsprl) then
    if (ispn.eq.1) then
      vl(:)=vl(:)+0.5d0*vqlss(:)
      vc(:)=vc(:)+0.5d0*vqcss(:)
    else
      vl(:)=vl(:)-0.5d0*vqlss(:)
      vc(:)=vc(:)-0.5d0*vqcss(:)
    end if
  end if
  call gengpvec(vl,vc,ngp(ispn),igpig(:,ispn),vgpl(:,:,ispn),vgpc(:,:,ispn))
  vl(:)=vpql(:)
  vc(:)=vpqc(:)
  if (spinsprl) then
    if (ispn.eq.1) then
      vl(:)=vl(:)+0.5d0*vqlss(:)
      vc(:)=vc(:)+0.5d0*vqcss(:)
    else
      vl(:)=vl(:)-0.5d0*vqlss(:)
      vc(:)=vc(:)-0.5d0*vqcss(:)
    end if
  end if
  call gengpvec(vl,vc,ngpq(ispn),igpqig(:,ispn),vgpql(:,:,ispn),vgpqc(:,:,ispn))
end do
! get the plane wave wavefunctions at p and p+q
call getwfpw(vkl(:,ikp),vgpl,wfpwp)
call getwfpw(vpql,vgpql,wfpwpq)
! compute the wavefunctions in real space
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispn,jspn,igp,igpq,ifg)
!$OMP DO
do ist=1,nstsv
  do ispn=1,nspinor
    if (spinsprl) then
      jspn=ispn
    else
      jspn=1
    end if
    wfirp(:,ispn,ist)=0.d0
    do igp=1,ngp(jspn)
      ifg=igfft(igpig(igp,jspn))
      wfirp(ifg,ispn,ist)=wfpwp(igp,ispn,ist)
    end do
    call zfftifc(3,ngrid,1,wfirp(:,ispn,ist))
    wfirpq(:,ispn,ist)=0.d0
    do igpq=1,ngpq(jspn)
      ifg=igfft(igpqig(igpq,jspn))
      wfirpq(ifg,ispn,ist)=wfpwpq(igpq,ispn,ist)
    end do
    call zfftifc(3,ngrid,1,wfirpq(:,ispn,ist))
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
t0=wkptnr/omega
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zw,zrhoir,zv,jst,eij) &
!$OMP PRIVATE(t1,iw,ig,jg,zt1,zt2)
!$OMP DO
do ist=1,nstsv
  allocate(zw(nwrpa))
  allocate(zrhoir(ngrtot))
  allocate(zv(ngrpa))
  do jst=1,nstsv
    eij=evalsv(ist,jkp)-(evalsv(jst,jkpq)+scissor)
! coefficient in response function formula for all RPA frequencies
! (note: this formula is unsuitable for systems without time-reversal symmetry)
    t1=t0*occsv(ist,jkp)*(1.d0-occsv(jst,jkpq)/occmax)
    do iw=1,nwrpa
      zw(iw)=t1*(1.d0/(eij+wrpa(iw))+1.d0/(eij-wrpa(iw)))
    end do
! compute the complex density in real space
    if (spinpol) then
      zrhoir(:)=conjg(wfirpq(:,1,jst))*wfirp(:,1,ist) &
               +conjg(wfirpq(:,2,jst))*wfirp(:,2,ist)
    else
      zrhoir(:)=conjg(wfirpq(:,1,jst))*wfirp(:,1,ist)
    end if
! Fourier transform density to G-space
    call zfftifc(3,ngrid,-1,zrhoir)
    do ig=1,ngrpa
      zv(ig)=zrhoir(igfft(ig))
    end do
! add to the response function (in general chi is not Hermitian)
!$OMP CRITICAL
    do ig=1,ngrpa
      zt1=conjg(zv(ig))
      do jg=1,ngrpa
        zt2=zt1*zv(jg)
        chi0(:,ig,jg)=chi0(:,ig,jg)+zt2*zw(:)
      end do
    end do
!$OMP END CRITICAL
  end do
  deallocate(zw,zrhoir,zv)
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(igpig,igpqig,vgpl,vgpql,vgpc,vgpqc)
deallocate(wfpwp,wfpwpq,wfirp,wfirpq)
return
end subroutine

