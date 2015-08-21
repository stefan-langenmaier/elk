
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlrad
! !INTERFACE:
subroutine hmlrad
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the radial Hamiltonian integrals of the APW and local-orbital
!   basis functions. In other words, for atom $\alpha$, it computes integrals of
!   the form
!   $$ h^{\alpha}_{qq';ll'l''m''}=\begin{cases}
!    \int_0^{R_i}u^{\alpha}_{q;l}(r)H u^{\alpha}_{q';l'}(r)r^2dr & l''=0 \\
!    \int_0^{R_i}u^{\alpha}_{q;l}(r)V^{\alpha}_{l''m''}(r)
!    u^{\alpha}_{q';l'}(r)r^2dr & l''>0 \end{cases}, $$
!   where $u^{\alpha}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; $H$ is the Hamiltonian of the radial Schr\"{o}dinger equation;
!   and $V^{\alpha}_{l''m''}$ is the muffin-tin Kohn-Sham potential. Similar
!   integrals are calculated for APW-local-orbital and
!   local-orbital-local-orbital contributions.
!
! !REVISION HISTORY:
!   Created December 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,nr,ir
integer l1,l2,l3,m2,lm2
integer ilo,jlo,io,jo
real(8) t1
! allocatable arrays
real(8), allocatable :: fr(:),gr(:)
! automatic arrays
real(8) r2(nrmtmax)
! begin loops over atoms and species
do is=1,nspecies
  nr=nrmt(is)
  r2(1:nr)=spr(1:nr,is)**2
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(fr,gr,ias) &
!$OMP PRIVATE(l1,l2,l3,io,jo,ir) &
!$OMP PRIVATE(lm2,m2,t1,ilo,jlo) SHARED(is)
!$OMP DO
  do ia=1,natoms(is)
    allocate(fr(nrmtmax),gr(nrmtmax))
    ias=idxas(ia,is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
    do l1=0,lmaxmat
      do io=1,apword(l1,is)
        do l3=0,lmaxmat
          do jo=1,apword(l3,is)
            if (l1.eq.l3) then
              do ir=1,nr
                fr(ir)=apwfr(ir,1,io,l1,ias)*apwfr(ir,2,jo,l3,ias)*r2(ir)
              end do
              call fderiv(-1,nr,spr(:,is),fr,gr)
              haa(1,jo,l3,io,l1,ias)=gr(nr)/y00
            else
              haa(1,jo,l3,io,l1,ias)=0.d0
            end if
            if (l1.ge.l3) then
              lm2=1
              do l2=1,lmaxvr
                do m2=-l2,l2
                  lm2=lm2+1
                  do ir=1,nr
                    t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*r2(ir)
                    fr(ir)=t1*vsmt(lm2,ir,ias)
                  end do
                  call fderiv(-1,nr,spr(:,is),fr,gr)
                  haa(lm2,jo,l3,io,l1,ias)=gr(nr)
                  haa(lm2,io,l1,jo,l3,ias)=gr(nr)
                end do
              end do
            end if
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
          if (l1.eq.l3) then
            do ir=1,nr
              fr(ir)=lofr(ir,1,ilo,ias)*apwfr(ir,2,io,l3,ias)*r2(ir)
            end do
            call fderiv(-1,nr,spr(:,is),fr,gr)
            hloa(1,io,l3,ilo,ias)=gr(nr)/y00
          else
            hloa(1,io,l3,ilo,ias)=0.d0
          end if
          lm2=1
          do l2=1,lmaxvr
            do m2=-l2,l2
              lm2=lm2+1
              do ir=1,nr
                t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*r2(ir)
                fr(ir)=t1*vsmt(lm2,ir,ias)
              end do
              call fderiv(-1,nr,spr(:,is),fr,gr)
              hloa(lm2,io,l3,ilo,ias)=gr(nr)
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
        if (l1.eq.l3) then
          do ir=1,nr
            fr(ir)=lofr(ir,1,ilo,ias)*lofr(ir,2,jlo,ias)*r2(ir)
          end do
          call fderiv(-1,nr,spr(:,is),fr,gr)
          hlolo(1,jlo,ilo,ias)=gr(nr)/y00
        else
          hlolo(1,jlo,ilo,ias)=0.d0
        end if
        lm2=1
        do l2=1,lmaxvr
          do m2=-l2,l2
            lm2=lm2+1
            do ir=1,nr
              t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*r2(ir)
              fr(ir)=t1*vsmt(lm2,ir,ias)
            end do
            call fderiv(-1,nr,spr(:,is),fr,gr)
            hlolo(lm2,jlo,ilo,ias)=gr(nr)
          end do
        end do
      end do
    end do
    deallocate(fr,gr)
! end loops over atoms and species
  end do
!$OMP END DO
!$OMP END PARALLEL
end do
return
end subroutine
!EOC

