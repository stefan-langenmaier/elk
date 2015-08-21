
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfpack(tpack,n,nr,nri,ld,zfmt,zfir,v)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(inout) :: n
integer, intent(in) :: nr(nspecies),nri(nspecies)
integer, intent(in) :: ld
complex(8), intent(inout) :: zfmt(lmmaxvr,ld,natmtot)
complex(8), intent(inout) :: zfir(ngtot)
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
        v(n+1)=dble(zfmt(lm,ir,ias))
        v(n+2)=aimag(zfmt(lm,ir,ias))
        n=n+2
      end do
      if (ir.eq.nri(is)) lmmax=lmmaxvr
    end do
  end do
  call dcopy(2*ngtot,zfir,1,v(n+1),1)
  n=n+2*ngtot
else
! unpack the function
  do ias=1,natmtot
    is=idxis(ias)
    lmmax=lmmaxinr
    do ir=1,nr(is)
      do lm=1,lmmax
        zfmt(lm,ir,ias)=cmplx(v(n+1),v(n+2),8)
        n=n+2
      end do
      if (ir.eq.nri(is)) lmmax=lmmaxvr
    end do
  end do
  call dcopy(2*ngtot,v(n+1),1,zfir,1)
  n=n+2*ngtot
end if
return
end subroutine

