
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlaa
! !INTERFACE:
subroutine hmlaa(ias,ngp,apwalm,h)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ias    : joint atom and species number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   h      : Hamiltonian matrix (inout,complex(*))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ias
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(inout) :: h(*)
! local variables
integer is,ld,io,jo
integer l1,l2,l3,m1,m2,m3
integer lm1,lm2,lm3
real(8) t1
complex(8) zt1,zsum
! automatic arrays
complex(8) zv(ngp)
ld=ngp+nlotot
is=idxis(ias)
lm1=0
do l1=0,lmaxmat
  do m1=-l1,l1
    lm1=lm1+1
    do io=1,apword(l1,is)
      zv(:)=0.d0
      lm3=0
      do l3=0,lmaxmat
        do m3=-l3,l3
          lm3=lm3+1
          if (lm1.ge.lm3) then
            do jo=1,apword(l3,is)
              zsum=0.d0
              do l2=0,lmaxvr
                if (mod(l1+l2+l3,2).eq.0) then
                  do m2=-l2,l2
                    lm2=idxlm(l2,m2)
                    zsum=zsum+gntyry(lm1,lm2,lm3)*haa(lm2,jo,l3,io,l1,ias)
                  end do
                end if
              end do
              if (lm1.eq.lm3) zsum=zsum*0.5d0
              if (abs(dble(zsum))+abs(aimag(zsum)).gt.1.d-14) then
                call zaxpy(ngp,zsum,apwalm(:,jo,lm3,ias),1,zv,1)
              end if
            end do
          end if
        end do
      end do
      call zher2b(ngp,zone,apwalm(:,io,lm1,ias),zv,ld,h)
    end do
  end do
end do
! kinetic surface contribution
t1=0.25d0*rmt(is)**2
lm1=0
do l1=0,lmaxmat
  do m1=-l1,l1
    lm1=lm1+1
    do io=1,apword(l1,is)
      do jo=1,apword(l1,is)
        zt1=t1*apwfr(nrmt(is),1,io,l1,ias)*apwdfr(jo,l1,ias)
        call zher2b(ngp,zt1,apwalm(:,io,lm1,ias),apwalm(:,jo,lm1,ias),ld,h)
      end do
    end do
  end do
end do
return
end subroutine
!EOC

