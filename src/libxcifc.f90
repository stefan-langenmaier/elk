
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module libxcifc

use xc_f90_lib_m

contains

!BOP
! !ROUTINE: xcifc_libxc
! !INTERFACE:
subroutine xcifc_libxc(xctype,n,rho,rhoup,rhodn,grho2,gup2,gdn2,gupdn,ex,ec, &
 vx,vc,vxup,vxdn,vcup,vcdn,dxdg2,dxdgu2,dxdgd2,dxdgud,dcdg2,dcdgu2,dcdgd2, &
 dcdgud)
! !INPUT/OUTPUT PARAMETERS:
!   xctype : type of exchange-correlation functional (in,integer(3))
!   n      : number of density points (in,integer)
!   rho    : spin-unpolarised charge density (in,real(n),optional)
!   rhoup  : spin-up charge density (in,real(n),optional)
!   rhodn  : spin-down charge density (in,real(n),optional)
!   grho2  : |grad rho|^2 (in,real(n),optional)
!   gup2   : |grad rhoup|^2 (in,real(n),optional)
!   gdn2   : |grad rhodn|^2 (in,real(n),optional)
!   gupdn  : (grad rhoup).(grad rhodn) (in,real(n),optional)
!   ex     : exchange energy density (out,real(n),optional)
!   ec     : correlation energy density (out,real(n),optional)
!   vx     : spin-unpolarised exchange potential (out,real(n),optional)
!   vc     : spin-unpolarised correlation potential (out,real(n),optional)
!   vxup   : spin-up exchange potential (out,real(n),optional)
!   vxdn   : spin-down exchange potential (out,real(n),optional)
!   vcup   : spin-up correlation potential (out,real(n),optional)
!   vcdn   : spin-down correlation potential (out,real(n),optional)
!   dxdg2  : de_x/d(|grad rho|^2) (out,real(n),optional)
!   dxdgu2 : de_x/d(|grad rhoup|^2) (out,real(n),optional)
!   dxdgd2 : de_x/d(|grad rhodn|^2) (out,real(n),optional)
!   dxdgud : de_x/d((grad rhoup).(grad rhodn)) (out,real(n),optional)
!   dcdg2  : de_c/d(|grad rho|^2) (out,real(n),optional)
!   dcdgu2 : de_c/d(|grad rhoup|^2) (out,real(n),optional)
!   dcdgd2 : de_c/d(|grad rhodn|^2) (out,real(n),optional)
!   dcdgud : de_c/d((grad rhoup).(grad rhodn)) (out,real(n),optional)
! !DESCRIPTION:
!   Interface to the ETSF {\tt libxc} exchange-correlation functional library:
!   \newline{\tt http://www.tddft.org/programs/octopus/wiki/index.php/Libxc}.
!   The second and third integers in {\tt xctype} define the exchange and
!   correlation functionals in {\tt libxc}, respectively.
!
! !REVISION HISTORY:
!   Created April 2009 (Tyrel McQueen)
!   Modified September 2009 (JKD and TMQ)
!   Updated for the libxc 1.0 interface, July 2010 (JKD)
!EOP
!BOC
implicit none
! mandatory arguments
integer, intent(in) :: xctype(3)
integer, intent(in) :: n
! optional arguments
real(8), optional, intent(in) :: rho(n)
real(8), optional, intent(in) :: rhoup(n)
real(8), optional, intent(in) :: rhodn(n)
real(8), optional, intent(in) :: grho2(n)
real(8), optional, intent(in) :: gup2(n)
real(8), optional, intent(in) :: gdn2(n)
real(8), optional, intent(in) :: gupdn(n)
real(8), optional, intent(out) :: ex(n)
real(8), optional, intent(out) :: ec(n)
real(8), optional, intent(out) :: vx(n)
real(8), optional, intent(out) :: vc(n)
real(8), optional, intent(out) :: vxup(n)
real(8), optional, intent(out) :: vxdn(n)
real(8), optional, intent(out) :: vcup(n)
real(8), optional, intent(out) :: vcdn(n)
real(8), optional, intent(out) :: dxdg2(n)
real(8), optional, intent(out) :: dxdgu2(n)
real(8), optional, intent(out) :: dxdgd2(n)
real(8), optional, intent(out) :: dxdgud(n)
real(8), optional, intent(out) :: dcdg2(n)
real(8), optional, intent(out) :: dcdgu2(n)
real(8), optional, intent(out) :: dcdgd2(n)
real(8), optional, intent(out) :: dcdgud(n)
! local variables
integer nspin,xcf,id,k
type(xc_f90_pointer_t) p,info
! allocatable arrays
real(8), allocatable :: r(:,:),sigma(:,:),v(:,:),vsigma(:,:)
if (present(rho)) then
  nspin=XC_UNPOLARIZED
else if (present(rhoup).and.present(rhodn)) then
  nspin=XC_POLARIZED
else
  write(*,*)
  write(*,'("Error(xcifc_libxc): missing arguments")')
  write(*,*)
  stop
end if
! loop over functional kinds (exchange or correlation)
do k=2,3
  id=xctype(k)
  if (id.gt.0) then
    xcf=xc_f90_family_from_id(id)
    select case(xcf)
    case(XC_FAMILY_LDA)
!-------------------------!
!     LDA functionals     !
!-------------------------!
      call xc_f90_func_init(p,info,id,nspin)
      if (k.eq.2) then
! exchange
        if (present(rho)) then
          call xc_f90_lda_exc_vxc(p,n,rho(1),ex(1),vx(1))
        else
          allocate(r(2,n),v(2,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          call xc_f90_lda_exc_vxc(p,n,r(1,1),ex(1),v(1,1))
          vxup(:)=v(1,:); vxdn(:)=v(2,:)
          deallocate(r,v)
        end if
      else
! correlation
        if (present(rho)) then
          call xc_f90_lda_exc_vxc(p,n,rho(1),ec(1),vc(1))
        else
          allocate(r(2,n),v(2,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          call xc_f90_lda_exc_vxc(p,n,r(1,1),ec(1),v(1,1))
          vcup(:)=v(1,:); vcdn=v(2,:)
          deallocate(r,v)
        end if
      end if
! destroy functional
      call xc_f90_func_end(p)
    case(XC_FAMILY_GGA)
!-------------------------!
!     GGA functionals     !
!-------------------------!
      call xc_f90_func_init(p,info,id,nspin)
      if (k.eq.2) then
! exchange
        if (present(rho)) then
          call xc_f90_gga_exc_vxc(p,n,rho(1),grho2(1),ex(1),vx(1),dxdg2(1))
        else
          allocate(r(2,n),sigma(3,n),v(2,n),vsigma(3,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          sigma(1,:)=gup2(:); sigma(2,:)=gupdn(:); sigma(3,:)=gdn2(:)
          call xc_f90_gga_exc_vxc(p,n,r(1,1),sigma(1,1),ex(1),v(1,1), &
           vsigma(1,1))
          vxup(:)=v(1,:); vxdn(:)=v(2,:)
          dxdgu2(:)=vsigma(1,:); dxdgud(:)=vsigma(2,:); dxdgd2(:)=vsigma(3,:)
          deallocate(r,sigma,v,vsigma)
        end if
      else
! correlation
        if (present(rho)) then
          call xc_f90_gga_exc_vxc(p,n,rho(1),grho2(1),ec(1),vc(1),dcdg2(1))
        else
          allocate(r(2,n),sigma(3,n),v(2,n),vsigma(3,n))
          r(1,:)=rhoup(:); r(2,:)=rhodn(:)
          sigma(1,:)=gup2(:); sigma(2,:)=gupdn(:); sigma(3,:)=gdn2(:)
          call xc_f90_gga_exc_vxc(p,n,r(1,1),sigma(1,1),ec(1),v(1,1), &
           vsigma(1,1))
          vcup(:)=v(1,:); vcdn(:)=v(2,:)
          dcdgu2(:)=vsigma(1,:); dcdgud(:)=vsigma(2,:); dcdgd2(:)=vsigma(3,:)
          deallocate(r,sigma,v,vsigma)
        end if
      end if
    case default
      write(*,*)
      write(*,'("Error(xcifc_libxc): unsupported libxc functional family : ",&
       &I8)') xcf
      write(*,*)
      stop
    end select
  else
! case when id=0
    if (k.eq.2) then
      ex(1:n)=0.d0
      if (present(rho)) then
        vx(1:n)=0.d0
      else
        vxup(1:n)=0.d0
        vxdn(1:n)=0.d0
      end if
    else
      ec(1:n)=0.d0
      if (present(rho)) then
        vc(1:n)=0.d0
      else
        vcup(1:n)=0.d0
        vcdn(1:n)=0.d0
      end if
    end if
  end if
end do
return
end subroutine

subroutine xcdata_libxc(xctype,xcdescr,xcspin,xcgrad)
implicit none
! arguments
integer, intent(in) :: xctype(3)
character(512), intent(out) :: xcdescr
integer, intent(out) :: xcspin
integer, intent(out) :: xcgrad
! local variables
integer xcf,id,k
character(256) name
type(xc_f90_pointer_t) p,info
! unknown spin polarisation
xcspin=-1
! no gradients by default
xcgrad=0
do k=2,3
  id=xctype(k)
  if (id.gt.0) then
    xcf=xc_f90_family_from_id(id)
    if ((xcf.ne.XC_FAMILY_LDA).and.(xcf.ne.XC_FAMILY_GGA)) then
      write(*,*)
      write(*,'("Error(xcdata_libxc): unsupported libxc functional family : ",&
       &I8)') xcf
      write(*,*)
      stop
    end if
! post-processed gradients required
    if (xcf.eq.XC_FAMILY_GGA) xcgrad=2
    call xc_f90_func_init(p,info,id,XC_UNPOLARIZED)
    call xc_f90_info_name(info,name)
    call xc_f90_func_end(p)
  else
    name='none'
  end if
  if (k.eq.2) then
    xcdescr='libxc; exchange: '//trim(name)
  else
    xcdescr=trim(xcdescr)//'; correlation: '//trim(name)
  end if
end do
xcdescr=trim(xcdescr)//' (see libxc for references)'
return
end subroutine
!EOC

end module

