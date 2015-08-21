
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gencore
! !INTERFACE:
subroutine gencore
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the core radial wavefunctions, eigenvalues and densities. The
!   radial Dirac equation is solved in the spherical part of the Kohn-Sham
!   potential to which the atomic potential has been appended for
!   $r>R_{\rm MT}$. In the case of spin-polarised calculations, and when
!   {\tt spincore} is set to {\tt .true.}, the Dirac equation is solved in the
!   spin-up and -down potentials created from the Kohn-Sham scalar potential and
!   magnetic field magnitude, with the occupancy divided equally between up and
!   down. The up and down densities determined in this way are added to both the
!   scalar density and the magnetisation in the routine {\tt rhocore}. Note
!   that this procedure is a simple, but inexact, approach to solving the radial
!   Dirac equation in a magnetic field.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Added polarised cores, November 2009 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,jas
integer ispn,ist,nr,ir
real(8) t1
! automatic arrays
logical done(natmmax)
real(8) vr(spnrmax),eval(spnstmax)
! allocatable arrays
real(8), allocatable :: br(:)
if (spincore) allocate(br(nrmtmax))
! loop over species and atoms
do is=1,nspecies
  nr=nrmt(is)
  done(:)=.false.
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
! Kohn-Sham magnetic field for spin-polarised core
    if (spincore) then
      if (ncmag) then
        do ir=1,nr
          br(ir)=sqrt(bxcmt(1,ir,ias,1)**2 &
                     +bxcmt(1,ir,ias,2)**2 &
                     +bxcmt(1,ir,ias,3)**2)*y00
        end do
      else
        do ir=1,nr
          br(ir)=abs(bxcmt(1,ir,ias,1))*y00
        end do
      end if
    end if
! loop over spin channels
    do ispn=1,nspncr
      if (frozencr) then
! use atomic potential for the frozen core approximation
        vr(1:nr)=spvr(1:nr,is)
      else
! else use the spherical part of the crystal Kohn-Sham potential
        vr(1:nr)=vsmt(1,1:nr,ias)*y00
      end if
! spin-up and -down potentials for polarised core
      if (spincore) then
        if (ispn.eq.1) then
          vr(1:nr)=vr(1:nr)+br(1:nr)
        else
          vr(1:nr)=vr(1:nr)-br(1:nr)
        end if
      end if
! append the Kohn-Sham potential from the atomic calculation for r > R_MT
      t1=vr(nr)-spvr(nr,is)
      do ir=nr+1,spnr(is)
        vr(ir)=spvr(ir,is)+t1
      end do
      rhocr(:,ias,ispn)=0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(t1,ir) SHARED(is)
!$OMP DO
      do ist=1,spnst(is)
        if (spcore(ist,is)) then
! solve the Dirac equation
          eval(ist)=evalcr(ist,ias)
          call rdirac(solsc,spn(ist,is),spl(ist,is),spk(ist,is),spnr(is), &
           spr(:,is),vr,eval(ist),rwfcr(:,1,ist,ias),rwfcr(:,2,ist,ias))
          if (spincore) then
! use the spin-averaged eigenvalue for the polarised core
            if (ispn.eq.1) then
              evalcr(ist,ias)=eval(ist)
            else
              evalcr(ist,ias)=0.5d0*(evalcr(ist,ias)+eval(ist))
            end if
            t1=0.5d0*occcr(ist,ias)
          else
            evalcr(ist,ias)=eval(ist)
            t1=occcr(ist,ias)
          end if
! add to the core density
!$OMP CRITICAL
          do ir=1,spnr(is)
            rhocr(ir,ias,ispn)=rhocr(ir,ias,ispn) &
             +t1*(rwfcr(ir,1,ist,ias)**2+rwfcr(ir,2,ist,ias)**2)
          end do
!$OMP END CRITICAL
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL
      do ir=1,spnr(is)
        rhocr(ir,ias,ispn)=rhocr(ir,ias,ispn)/(fourpi*spr(ir,is)**2)
      end do
! end loop over spin channels
    end do
    done(ia)=.true.
! copy to equivalent atoms
    do ja=1,natoms(is)
      if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
        jas=idxas(ja,is)
        do ist=1,spnst(is)
          if (spcore(ist,is)) then
            evalcr(ist,jas)=evalcr(ist,ias)
            rwfcr(1:spnr(is),:,ist,jas)=rwfcr(1:spnr(is),:,ist,ias)
          end if
        end do
        rhocr(1:spnr(is),jas,:)=rhocr(1:spnr(is),ias,:)
        done(ja)=.true.
      end if
    end do
! end loop over species and atoms
  end do
end do
if (spincore) deallocate(br)
return
end subroutine
!EOC

