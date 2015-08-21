
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dyntask(fnum)
use modmain
use modphonon
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
logical exist
do ipph=1,3
  do isph=1,nspecies
    do iaph=1,natoms(isph)
      do iqph=1,nqpt
! construct the phonon file extension
        call phfext(iqph,isph,iaph,ipph,filext)
! determine if the DYN file with this extension exists
        inquire(file='DYN'//trim(filext),exist=exist)
        if (.not.exist) then
          open(fnum,file='DYN'//trim(filext),action='WRITE',form='FORMATTED')
          return
        end if
      end do
    end do
  end do
end do
filext='.OUT'
iqph=0; isph=0; iaph=0; ipph=0
write(*,*)
write(*,'("Info(dyntask): nothing more to do")')
return
end subroutine

