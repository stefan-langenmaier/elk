
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: mossbauer
! !INTERFACE:
subroutine mossbauer
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the contact charge density and contact magnetic hyperfine field for
!   each atom and outputs the data to the file {\tt MOSSBAUER.OUT}. The nuclear
!   radius used for the contact density is approximated by the empirical formula
!   $R_{\rm N}=1.25 Z^{1/3}$ fm, where $Z$ is the atomic number. The Thomson
!   radius, $R_{\rm T}=Z/c^2$, is used for determining the contact moment $m_c$,
!   the relation of which to field strength is via Fermi's formula
!   $$ B_c=\frac{8\pi}{3}\mu_B m_c, $$
!   where the orbital and dipolar contributions are neglected. See
!   S. Bl\"{u}gel, H. Akai, R. Zeller, and P. H. Dederichs, {\it Phys. Rev. B}
!   {\bf 35}, 3271 (1987).
!
! !REVISION HISTORY:
!   Created May 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias
integer ir,nrn,nrt
! nuclear radius constant in Bohr
real(8), parameter :: r0=1.25d-15/au_si
real(8) rn,rt,vn,vt
real(8) rho0,mc,bc,t1
! allocatable arrays
real(8), allocatable :: fr(:)
real(8), allocatable :: gr(:)
real(8), allocatable :: cf(:,:)
! initialise universal variables
call init0
! read density and potentials from file
call readstate
! allocate local arrays
allocate(fr(nrmtmax))
allocate(gr(nrmtmax))
allocate(cf(4,nrmtmax))
open(50,file='MOSSBAUER.OUT',action='WRITE',form='FORMATTED')
do is=1,nspecies
! approximate nuclear radius and volume
  rn=r0*abs(spzn(is))**(1.d0/3.d0)
  do ir=1,nrmt(is)-1
    if (spr(ir,is).gt.rn) goto 10
  end do
10 continue
  nrn=ir
  rn=spr(nrn,is)
  vn=(4.d0/3.d0)*pi*rn**3
! Thomson radius and volume
  rt=abs(spzn(is))/solsc**2
  do ir=1,nrmt(is)-1
    if (spr(ir,is).gt.rt) goto 20
  end do
20 continue
  nrt=ir
  rt=spr(nrt,is)
  vt=(4.d0/3.d0)*pi*rt**3
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!--------------------------------!
!     contact charge density     !
!--------------------------------!
    do ir=1,nrn
      fr(ir)=(fourpi*spr(ir,is)**2)*y00*rhomt(1,ir,ias)
    end do
    call fderiv(-1,nrn,spr(:,is),fr,gr,cf)
    rho0=gr(nrn)/vn
    write(50,*)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    write(50,*)
    write(50,'(" approximate nuclear radius : ",G18.10)') rn
    write(50,'(" number of mesh points to nuclear radius : ",I6)') nrn
    write(50,'(" average contact charge density : ",G18.10)') rho0
!------------------------------------------!
!     contact magnetic hyperfine field     !
!------------------------------------------!
    if (spinpol) then
      do ir=1,nrt
        if (ncmag) then
! non-collinear
          t1=sqrt(magmt(1,ir,ias,1)**2 &
                 +magmt(1,ir,ias,2)**2 &
                 +magmt(1,ir,ias,3)**2)
        else
! collinear
          t1=magmt(1,ir,ias,1)
        end if
        fr(ir)=(fourpi*spr(ir,is)**2)*y00*t1
      end do
      call fderiv(-1,nrt,spr(:,is),fr,gr,cf)
      mc=gr(nrt)/vt
      write(50,*)
      write(50,'(" Thomson radius : ",G18.10)') rt
      write(50,'(" number of mesh points to Thomson radius : ",I6)') nrt
      write(50,'(" contact magnetic moment (mu_B) : ",G18.10)') mc
      bc=(8.d0*pi/3.d0)*mc/(2.d0*solsc)
      write(50,'(" contact hyperfine field : ",G18.10)') bc
      write(50,'("  tesla                  : ",G18.10)') bc*b_si/solsc
    end if
  end do
end do
close(50)
write(*,*)
write(*,'("Info(mossbauer):")')
write(*,'(" Mossbauer parameters written to MOSSBAUER.OUT")')
deallocate(fr,gr,cf)
return
end subroutine
!EOC
