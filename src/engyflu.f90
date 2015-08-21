
! Copyright (C) 2008 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: engyflu
! !INTERFACE:
subroutine engyflu
! !USES:
use modmain
use modldapu
! !DESCRIPTION:
!   Calculates the energies of radial functions to be used to calculate the
!   Slater integrals. By convention those energies are chosen to be the ones at
!   the center of the band.
!
! !REVISION HISTORY:
!   Created April 2008 from linengy (Francesco Cricchio)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,jas
integer l,io
logical fnd
! automatic arrays
logical done(natmmax)
real(8) vr(nrmtmax)
! begin loops over atoms and species
do is=1,nspecies
  l=llu(is)
  io=1
  if (l.lt.0) goto 10
  done(:)=.false.
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
      vr(1:nrmt(is))=veffmt(1,1:nrmt(is),ias)*y00
! find the center of the band starting from -0.5 htr
! in my experience this value is low enough in order to not miss it
! and jump by mistake, for a given l, at higher n quantum number
      flue(io,l,ias)=-0.5
      call findband(solsc,l,0,nprad,nrmt(is),spr(1,is),vr,deband,epsband, &
       flue(io,l,ias),fnd)
      if (.not.fnd) then
        write(*,*)
        write(*,'("Warning(engyflu): energy not found")')
        write(*,'(" for species ",I4)') is
        write(*,'(" atom ",I4)') ia
        write(*,'(" angular momentum ",I4)') l
        write(*,'(" and s.c. loop ",I5)') iscl
      end if
      done(ia)=.true.
! copy to equivalent atoms
      do ja=1,natoms(is)
        if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
          jas=idxas(ja,is)
          io=1
          flue(io,l,jas)=flue(io,l,ias)
          done(ja)=.true.
        end if
      end do
    end if
! end loops over atoms and species
  end do
10 continue
end do
return
end subroutine
!EOC

