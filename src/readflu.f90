
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readflu
use modmain
use modldapu
implicit none
! local variables
integer is,ia,l,k
flu(:,:)=0.d0
! read Slater integrals from FLU.OUT
open(50,file='FLU'//trim(filext),action='READ',form='FORMATTED')
do is=1,nspecies
  l=llu(is)
  if (l.ge.0) then
    do ia=1,natoms(is)
      read(50,*)
      read(50,*)
      do k=0,2*l,2
        read(50,*) flu(k,is)
      end do
      if (inptypelu.ge.4) read(50,*)
      read(50,*)
      read(50,*)
    end do
  end if
end do
close(50)
return
end subroutine

