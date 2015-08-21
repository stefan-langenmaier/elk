
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine mixander(iscl,beta,n,nu,mu,f,g,d)
implicit none
! arguments
integer, intent(in) :: iscl
real(8), intent(in) :: beta
integer, intent(in) :: n
real(8), intent(inout) :: nu(n)
real(8), intent(inout) :: mu(n,3)
real(8), intent(inout) :: f(n,3)
real(8), intent(inout) :: g(n,2)
real(8), intent(out) :: d
! local variables
integer j0,j1,j2,k
real(8) d11,d12,d22,r1,r2,t1,t2
if (iscl.le.1) then
  mu(:,1)=nu(:)
  f(:,1)=0.d0
  d=1.d0
  return
end if
j0=mod(iscl-1,3)+1
j1=mod(iscl-2,3)+1
j2=mod(iscl-3,3)+1
d=0.d0
do k=1,n
  f(k,j0)=nu(k)-mu(k,j1)
  d=d+f(k,j0)**2
end do
d=sqrt(d/dble(n))
if (iscl.eq.2) then
  nu(:)=beta*nu(:)+(1.d0-beta)*mu(:,j1)
  mu(:,j0)=nu(:)
  return
end if
g(:,1)=f(:,j1)-f(:,j0)
d11=dot_product(g(:,1),g(:,1))
r1=dot_product(f(:,j0),g(:,1))
if (iscl.eq.3) then
  t1=-r1/d11
  nu(:)=(1.d0-t1)*mu(:,j1)+t1*mu(:,j2)+beta*(f(:,j0)+t1*g(:,1))
else
  g(:,2)=f(:,j2)-f(:,j0)
  d22=dot_product(g(:,2),g(:,2))
  d12=dot_product(g(:,1),g(:,2))
  r2=dot_product(f(:,j0),g(:,2))
  t1=d11*d22
  t2=t1-d12*d12
  if (abs(t2).gt.abs(t1)*1.d-8) then
    t1=(-r1*d22+r2*d12)/t2
    t2=(r1*d12-r2*d11)/t2
  else
    t1=0.d0
    t2=0.d0
  end if
  nu(:)=(1.d0-t1-t2)*mu(:,j1)+t1*mu(:,j2)+t2*mu(:,j0)+beta*(f(:,j0)+t1*g(:,1) &
   +t2*g(:,2))
end if
mu(:,j0)=nu(:)
return
end subroutine

