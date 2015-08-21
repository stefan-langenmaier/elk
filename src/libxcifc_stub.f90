
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Stub routines for libxc. See Elk manual for libxc installation instructions.

module libxcifc

contains

subroutine xcifc_libxc(xctype,n,c_tb09,rho,rhoup,rhodn,g2rho,g2up,g2dn,grho2, &
 gup2,gdn2,gupdn,tau,ex,ec,vx,vc,vxup,vxdn,vcup,vcdn,dxdg2,dxdgu2,dxdgd2, &
 dxdgud,dcdg2,dcdgu2,dcdgd2,dcdgud)
implicit none
! mandatory arguments
integer xctype(3),n
! optional arguments
real(8), optional :: c_tb09
real(8), optional :: rho(n)
real(8), optional :: rhoup(n)
real(8), optional :: rhodn(n)
real(8), optional :: g2rho(n)
real(8), optional :: g2up(n)
real(8), optional :: g2dn(n)
real(8), optional :: grho2(n)
real(8), optional :: gup2(n)
real(8), optional :: gdn2(n)
real(8), optional :: gupdn(n)
real(8), optional :: tau(n)
real(8), optional :: ex(n)
real(8), optional :: ec(n)
real(8), optional :: vx(n)
real(8), optional :: vc(n)
real(8), optional :: vxup(n)
real(8), optional :: vxdn(n)
real(8), optional :: vcup(n)
real(8), optional :: vcdn(n)
real(8), optional :: dxdg2(n)
real(8), optional :: dxdgu2(n)
real(8), optional :: dxdgd2(n)
real(8), optional :: dxdgud(n)
real(8), optional :: dcdg2(n)
real(8), optional :: dcdgu2(n)
real(8), optional :: dcdgd2(n)
real(8), optional :: dcdgud(n)
write(*,*)
write(*,'("Error(libxcifc): libxc not or improperly installed")')
write(*,*)
stop
end subroutine

subroutine xcdata_libxc(xctype,xcdescr,xcspin,xcgrad)
implicit none
! arguments
integer :: xctype(3)
character(512) :: xcdescr
integer :: xcspin
integer :: xcgrad
write(*,*)
write(*,'("Error(libxcifc):  libxc not or improperly installed")')
write(*,*)
stop
end subroutine
!EOC

end module

