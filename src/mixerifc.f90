
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mixerifc(mtype,n,v,dv,nwork,work)
use modmain
implicit none
! arguments
integer, intent(in) :: mtype
integer, intent(in) :: n
real(8), intent(inout) :: v(n)
real(8), intent(out) :: dv
integer, intent(inout) :: nwork
real(8), intent(inout) :: work(*)
! local variables
! maximum subspace dimension for the Pulay mixer
integer, parameter :: maxsd=3
select case(mtype)
case(1)
! adaptive linear mixing
  if (nwork.le.0) then
    nwork=3*n
    return
  end if
  call mixadapt(iscl,beta0,betamax,n,v,work,work(n+1),work(2*n+1),dv)
case(2)
! Pulay mixing
  if (nwork.le.0) then
    nwork=2*maxsd*n
    return
  end if
  call mixpulay(iscl,n,maxsd,v,work,work(n*maxsd+1),dv)
case(3)
! Anderson mixing
  if (nwork.le.0) then
    nwork=8*n
    return
  end if
  call mixander(iscl,beta0,n,v,work,work(3*n+1),work(6*n+1),dv)
case default
  write(*,*)
  write(*,'("Error(mixerifc): mtype not defined : ",I8)') mtype
  write(*,*)
  stop
end select
return
end subroutine

subroutine getmixdata(mtype,mixdescr)
implicit none
! arguments
integer, intent(in) :: mtype
character(256), intent(out) :: mixdescr
select case(mtype)
case(1)
  mixdescr='Adaptive linear mixing'
case(2)
  mixdescr='Pulay mixing, Chem. Phys. Lett. 73, 393 (1980)'
case(3)
  mixdescr='Anderson mixing, J. Assoc. Comput. Mach. 12, 547 (1964)'
case default
  write(*,*)
  write(*,'("Error(getmixdata): mixtype not defined : ",I8)') mtype
  write(*,*)
  stop
end select
return
end subroutine

