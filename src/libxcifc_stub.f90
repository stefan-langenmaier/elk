
! Copyright (C) 2009 T. McQueen and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Stub routines for libxc. See Elk manual for libxc installation instructions.

module libxcifc

contains

subroutine xcifc_libxc(xctype,n,rho,rhoup,rhodn,grho2,gup2,gdn2,gupdn,ex,ec, &
 vx,vc,vxup,vxdn,vcup,vcdn,dxdg2,dxdgu2,dxdgd2,dxdgud,dcdg2,dcdgu2,dcdgd2, &
 dcdgud)
implicit none
integer xctype(3),n
real(8), optional :: rho(*),rhoup(*),rhodn(*)
real(8), optional :: grho2(*),gup2(*),gdn2(*),gupdn(*)
real(8), optional :: ex(*),ec(*),vx(*),vc(*)
real(8), optional :: vxup(*),vxdn(*),vcup(*),vcdn(*)
real(8), optional :: dxdg2(*),dxdgu2(*),dxdgd2(*),dxdgud(*)
real(8), optional :: dcdg2(*),dcdgu2(*),dcdgd2(*),dcdgud(*)
write(*,*)
write(*,'("Error(libxcifc): libxc not or improperly installed")')
write(*,*)
stop
end subroutine

subroutine xcdata_libxc(xctype,xcdescr,xcspin,xcgrad)
implicit none
integer xctype(3),xcspin,xcgrad
character(512) :: xcdescr
write(*,*)
write(*,'("Error(libxcifc):  libxc not or improperly installed")')
write(*,*)
stop
end subroutine
!EOC

end module

