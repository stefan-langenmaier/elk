
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfpack(tpack,n,nr,ld,rfmt,rfir,nu)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(inout) :: n
integer, intent(in) :: nr(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rfmt(lmmaxvr,ld,natmtot)
real(8), intent(inout) :: rfir(ngtot)
real(8), intent(out) :: nu(*)
! local variables
integer is,ias,ir,lm
if (tpack) then
! pack the function
  do ias=1,natmtot
    is=idxis(ias)
    do ir=1,nr(is)
      do lm=1,lmmaxvr
        n=n+1
        nu(n)=rfmt(lm,ir,ias)
      end do
    end do
  end do
  do ir=1,ngtot
    n=n+1
    nu(n)=rfir(ir)
  end do
else
! unpack the function
  do ias=1,natmtot
    is=idxis(ias)
    do ir=1,nr(is)
      do lm=1,lmmaxvr
        n=n+1
        rfmt(lm,ir,ias)=nu(n)
      end do
    end do
  end do
  do ir=1,ngtot
    n=n+1
    rfir(ir)=nu(n)
  end do
end if
return
end subroutine

