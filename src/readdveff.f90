
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readdveff(iq,is,ia,ip,dveffmt,dveffir)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ip
complex(8), intent(out) :: dveffmt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(out) :: dveffir(ngrtot)
! local variables
integer js,iostat
integer version_(3),nspecies_,lmmaxvr_
integer natoms_,nrcmt_,ngrid_(3)
character(256) fext,fname
call phfext(iq,is,ia,ip,fext)
fname='DVEFF'//trim(fext)
open(50,file=trim(fname),action='READ',form='UNFORMATTED',status='OLD', &
 iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readdveff): error opening ",A)') trim(fname)
  write(*,*)
  stop
end if
read(50) version_
if ((version(1).ne.version_(1)).or.(version(2).ne.version_(2)) &
 .or.(version(3).ne.version_(3))) then
  write(*,*)
  write(*,'("Warning(readdveff): different versions")')
  write(*,'(" current : ",I3.3,".",I3.3,".",I3.3)') version
  write(*,'(" file    : ",I3.3,".",I3.3,".",I3.3)') version_
  write(*,'(" in file ",A)') trim(fname)
end if
read(50) nspecies_
if (nspecies.ne.nspecies_) then
  write(*,*)
  write(*,'("Error(readdveff): differing nspecies")')
  write(*,'(" current : ",I4)') nspecies
  write(*,'(" file    : ",I4)') nspecies_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
read(50) lmmaxvr_
if (lmmaxvr.ne.lmmaxvr_) then
  write(*,*)
  write(*,'("Error(readdveff): differing lmmaxvr")')
  write(*,'(" current : ",I4)') lmmaxvr
  write(*,'(" file    : ",I4)') lmmaxvr_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
do js=1,nspecies
  read(50) natoms_
  if (natoms(js).ne.natoms_) then
    write(*,*)
    write(*,'("Error(readdveff): differing natoms for species ",I4)') js
    write(*,'(" current : ",I4)') natoms(js)
    write(*,'(" file    : ",I4)') natoms_
    write(*,'(" in file ",A)') trim(fname)
    write(*,*)
    stop
  end if
  read(50) nrcmt_
  if (nrcmt(js).ne.nrcmt_) then
    write(*,*)
    write(*,'("Error(readdveff): differing nrcmt for species ",I4)') js
    write(*,'(" current : ",I6)') nrcmt(js)
    write(*,'(" file    : ",I6)') nrcmt_
    write(*,'(" in file ",A)') trim(fname)
    write(*,*)
    stop
  end if
end do
read(50) ngrid_
if ((ngrid(1).ne.ngrid_(1)).or.(ngrid(2).ne.ngrid_(2)).or. &
 (ngrid(3).ne.ngrid_(3))) then
  write(*,*)
  write(*,'("Error(readdveff): differing ngrid")')
  write(*,'(" current : ",3I6)') ngrid
  write(*,'(" file    : ",3I6)') ngrid_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
read(50) dveffmt,dveffir
close(50)
return
end subroutine

