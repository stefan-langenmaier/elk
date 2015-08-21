
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdwfpw(vpl,dwfpw)
use modmain
use modphonon
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(out) :: dwfpw(ngkmax,nspinor,nstsv)
! local variables
integer isym,ikp,ikpq
integer ngp(nspnfv),ngpq(nspnfv)
integer igp,igpq,ifg,ir
integer ist,jst,ispn,jspn
real(8) vpc(3),vpql(3),vpqc(3)
real(8) vl(3),vc(3),eij
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
complex(8), allocatable :: wfir(:,:)
complex(8), allocatable :: zv(:)
! external functions
complex(8) zdotc
external zdotc
allocate(igpig(ngkmax,nspnfv))
allocate(igpqig(ngkmax,nspnfv))
allocate(vgpl(3,ngkmax,nspnfv))
allocate(vgpql(3,ngkmax,nspnfv))
allocate(vgpc(3,ngkmax,nspnfv))
allocate(vgpqc(3,ngkmax,nspnfv))
allocate(wfpwp(ngkmax,nspinor,nstsv))
allocate(wfpwpq(ngkmax,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor))
allocate(zv(ngkmax))
! p vector in Cartesian coordinates
call r3mv(bvec,vpl,vpc)
! p+q vector in lattice and Cartesian coordinates
vpql(:)=vpl(:)+vql(:,iqph)
vpqc(:)=vpc(:)+vqc(:,iqph)
! equivalent reduced k-points for p and p+q
call findkpt(vpl,isym,ikp)
call findkpt(vpql,isym,ikpq)
! generate the G+p and G+p+q vectors
do ispn=1,nspnfv
  vl(:)=vpl(:)
  vc(:)=vpc(:)
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
call getwfpw(vpl,vgpl,wfpwp)
call getwfpw(vpql,vgpql,wfpwpq)
! loop over states of k-point p
do ist=1,nstsv
  do ispn=1,nspinor
    if (spinsprl) then
      jspn=ispn
    else
      jspn=1
    end if
    wfir(:,ispn)=0.d0
    do igp=1,ngp(jspn)
      ifg=igfft(igpig(igp,jspn))
      wfir(ifg,ispn)=wfpwp(igp,ispn,ist)
    end do
! Fourier transform to real space
    call zfftifc(3,ngrid,1,wfir(:,ispn))
  end do
! multiply with the lattice-periodic part of the perturbing potential
! (i.e. without the exp(iq.r) phase)
  if (spinpol) then
! spin-polarised case
    if (ncmag) then
! non-collinear
      do ir=1,ngrtot
        zt1=cmplx(-aimag(dbxcpw(ir,2)),dble(dbxcpw(ir,2)),8)
        zt2=(dbxcpw(ir,1)+zt1)*wfir(ir,1) &
         +(dveffpw(ir)-dbxcpw(ir,3))*wfir(ir,2)
        wfir(ir,1)=(dveffpw(ir)+dbxcpw(ir,3))*wfir(ir,1) &
         +(dbxcpw(ir,1)-zt1)*wfir(ir,2)
        wfir(ir,2)=zt2
      end do
    else
! collinear
      do ir=1,ngrtot
        wfir(ir,1)=(dveffpw(ir)+dbxcpw(ir,3))*wfir(ir,1)
        wfir(ir,2)=(dveffpw(ir)-dbxcpw(ir,3))*wfir(ir,2)
      end do
    end if
  else
! spin-unpolarised case
    wfir(:,1)=wfir(:,1)*dveffpw(:)
  end if
  do ispn=1,nspinor
    if (spinsprl) then
      jspn=ispn
    else
      jspn=1
    end if
! Fourier transform to G-space
    call zfftifc(3,ngrid,-1,wfir(:,ispn))
    do igpq=1,ngpq(jspn)
      ifg=igfft(igpqig(igpq,jspn))
      zv(igpq)=wfir(ifg,ispn)
    end do
    dwfpw(:,ispn,ist)=0.d0
    do jst=1,nstsv
      eij=evalsv(ist,ikp)-evalsv(jst,ikpq)
      if (abs(eij).gt.1.d-8) then
        zt1=zdotc(ngpq(jspn),wfpwpq(:,ispn,jst),1,zv,1)
        call zaxpy(ngpq(jspn),zt1,wfpwpq(:,ispn,jst),1,dwfpw(:,ispn,ist),1)
      end if
    end do
  end do
end do
deallocate(igpig,igpqig)
deallocate(vgpl,vgpql,vgpc,vgpqc)
deallocate(wfpwp,wfpwpq,wfir,zv)
return
end subroutine

