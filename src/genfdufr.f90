
! Copyright (C) 2008  F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genfdufr
! !INTERFACE:
subroutine genfdufr
! !USES:
use modmain
use moddftu
! !DESCRIPTION:
!   Generates the radial functions used to calculate the Slater integrals
!   through a Yukawa potential.
!
! !REVISION HISTORY:
!   Created April 2008 from genapwfr (Francesco Cricchio)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,i
integer l,nn,io,nr,ir
real(8) t1
! automatic arrays
real(8) vr(nrmtmax),fr(nrmtmax),gr(nrmtmax)
real(8) p0(nrmtmax,apwordmax),p1(nrmtmax),p1s(apwordmax)
real(8) q0(nrmtmax,apwordmax),q1(nrmtmax,apwordmax)
real(8) hp0(nrmtmax)
do i=1,ndftu
  is=idftu(1,i)
  l=idftu(2,i)
  io=1
  nr=nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    vr(1:nr)=vsmt(1,1:nr,ias)*y00
! integrate the radial Schrodinger equation
    call rschroddme(solsc,0,l,0,fdue(io,l,ias),nr,spr(1,is),vr,nn,p0(1,io),p1, &
     q0(1,io),q1(1,io))
! normalise radial functions
    do ir=1,nr
      fr(ir)=p0(ir,io)**2
    end do
    call fderiv(-1,nr,spr(1,is),fr,gr)
    t1=1.d0/sqrt(abs(gr(nr)))
    p0(1:nr,io)=t1*p0(1:nr,io)
    p1s(io)=t1*p1(nr)
    q0(1:nr,io)=t1*q0(1:nr,io)
    q1(1:nr,io)=t1*q1(1:nr,io)
! normalise radial functions
    do ir=1,nr
      fr(ir)=p0(ir,io)**2
    end do
    call fderiv(-1,nr,spr(1,is),fr,gr)
    t1=abs(gr(nr))
    if (t1.lt.1.d-20) then
      write(*,*)
      write(*,'("Error(genfdufr): degenerate APW radial functions")')
      write(*,'(" for species ",I4)') is
      write(*,'(" atom ",I4)') ia
      write(*,'(" angular momentum ",I4)') l
      write(*,'(" and order ",I4)') io
      write(*,*)
      stop
    end if
    t1=1.d0/sqrt(t1)
    p0(1:nr,io)=t1*p0(1:nr,io)
    p1s(io)=t1*p1s(io)
    q0(1:nr,io)=t1*q0(1:nr,io)
    q1(1:nr,io)=t1*q1(1:nr,io)
! apply the Hamiltonian
    call rschrodapp(solsc,l,nr,spr(1,is),vr,p0(1,io),q0(1,io),q1(1,io),hp0)
! divide by r and store in global array
    do ir=1,nr
      t1=1.d0/spr(ir,is)
      fdufr(ir,1,io,l,ias)=t1*p0(ir,io)
      fdufr(ir,2,io,l,ias)=t1*hp0(ir)
    end do
  end do
end do
return
end subroutine
!EOC

