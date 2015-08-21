
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dhmlrad
use modmain
use modphonon
implicit none
! local variables
integer is,ia,ias,nr,ir
integer l1,l2,l3,m2,lm2
integer ilo,jlo,io,jo
real(8) t1,t2
! allocatable arrays
real(8), allocatable :: fr1(:),fr2(:),gr(:)
! automatic arrays
real(8) r2(nrmtmax)
! begin loops over atoms and species
do is=1,nspecies
  nr=nrmt(is)
  r2(1:nr)=spr(1:nr,is)**2
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(fr1,fr2,gr,ias) &
!$OMP PRIVATE(l1,l2,l3,io,jo,ir) &
!$OMP PRIVATE(lm2,m2,t1,t2,ilo,jlo) SHARED(is)
!$OMP DO
  do ia=1,natoms(is)
    allocate(fr1(nrmtmax),fr2(nrmtmax),gr(nrmtmax))
    ias=idxas(ia,is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
    do l1=0,lmaxmat
      do io=1,apword(l1,is)
        do l3=0,lmaxmat
          do jo=1,apword(l3,is)
            lm2=0
            do l2=0,lmaxvr
              do m2=-l2,l2
                lm2=lm2+1
                do ir=1,nr
                  t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*r2(ir)
                  fr1(ir)=t1*dble(dvsmt(lm2,ir,ias))
                  fr2(ir)=t1*aimag(dvsmt(lm2,ir,ias))
                end do
                call fderiv(-1,nr,spr(:,is),fr1,gr)
                t1=gr(nr)
                call fderiv(-1,nr,spr(:,is),fr2,gr)
                t2=gr(nr)
                dhaa(lm2,jo,l3,io,l1,ias)=cmplx(t1,t2,8)
              end do
            end do
          end do
        end do
      end do
    end do
!-------------------------------------!
!     local-orbital-APW integrals     !
!-------------------------------------!
    do ilo=1,nlorb(is)
      l1=lorbl(ilo,is)
      do l3=0,lmaxmat
        do io=1,apword(l3,is)
          lm2=0
          do l2=0,lmaxvr
            do m2=-l2,l2
              lm2=lm2+1
              do ir=1,nr
                t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*r2(ir)
                fr1(ir)=t1*dble(dvsmt(lm2,ir,ias))
                fr2(ir)=t1*aimag(dvsmt(lm2,ir,ias))
              end do
              call fderiv(-1,nr,spr(:,is),fr1,gr)
              t1=gr(nr)
              call fderiv(-1,nr,spr(:,is),fr2,gr)
              t2=gr(nr)
              dhloa(lm2,io,l3,ilo,ias)=cmplx(t1,t2,8)
            end do
          end do
        end do
      end do
    end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
    do ilo=1,nlorb(is)
      l1=lorbl(ilo,is)
      do jlo=1,nlorb(is)
        l3=lorbl(jlo,is)
        lm2=0
        do l2=0,lmaxvr
          do m2=-l2,l2
            lm2=lm2+1
            do ir=1,nr
              t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*r2(ir)
              fr1(ir)=t1*dble(dvsmt(lm2,ir,ias))
              fr2(ir)=t1*aimag(dvsmt(lm2,ir,ias))
            end do
            call fderiv(-1,nr,spr(:,is),fr1,gr)
            t1=gr(nr)
            call fderiv(-1,nr,spr(:,is),fr2,gr)
            t2=gr(nr)
            dhlolo(lm2,jlo,ilo,ias)=cmplx(t1,t2,8)
          end do
        end do
      end do
    end do
    deallocate(fr1,fr2,gr)
! end loops over atoms and species
  end do
!$OMP END DO
!$OMP END PARALLEL
end do
return
end subroutine

