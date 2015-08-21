
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfpack(tpack,n,nr,nri,ld,rfmt,rfir,v)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(inout) :: n
integer, intent(in) :: nr(nspecies),nri(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rfmt(lmmaxvr,ld,natmtot)
real(8), intent(inout) :: rfir(ngtot)
real(8), intent(out) :: v(*)
! local variables
integer is,ias,ir,lmmax,lm
if (tpack) then
! pack the function
  do ias=1,natmtot
    is=idxis(ias)
    lmmax=lmmaxinr
    do ir=1,nr(is)
      do lm=1,lmmax
        n=n+1
        v(n)=rfmt(lm,ir,ias)
      end do
      if (ir.eq.nri(is)) lmmax=lmmaxvr
    end do
  end do
  call dcopy(ngtot,rfir,1,v(n+1),1)
  n=n+ngtot
else
! unpack the function
  do ias=1,natmtot
    is=idxis(ias)
    lmmax=lmmaxinr
    do ir=1,nr(is)
      do lm=1,lmmax
        n=n+1
        rfmt(lm,ir,ias)=v(n)
      end do
      if (ir.eq.nri(is)) lmmax=lmmaxvr
    end do
  end do
  call dcopy(ngtot,v(n+1),1,rfir,1)
  n=n+ngtot
end if
return
end subroutine

