
! Copyright (C) 2009 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine readalphalu
use modmain
use modldapu
implicit none
! local variables
integer is,ia,ias,l
integer is_,ia_,l_
if (.not.readalu) return
open(50,file='ALPHALU'//trim(filext),action='READ',form='FORMATTED')
do is=1,nspecies
  l=llu(is)
  if (l.ge.0) then
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      read(50,*)
      read(50,*) is_,ia_,l_
      if ((is.ne.is_).or.(ia.ne.ia_).or.(l.ne.l_)) then
        write(*,*)
        write(*,'("Error(readalu): differing is, ia or l")')
        write(*,'(" current     : ",3I4)') is,ia,l
        write(*,'(" ALPHALU.OUT : ",3I4)') is_,ia_,l_
        write(*,*)
        stop
      end if
      read(50,*) alphalu(ias)
    end do
  end if
end do
close(50)
return
end subroutine

