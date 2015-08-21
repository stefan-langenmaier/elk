
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine energykncr
use modmain
implicit none
integer is,ia,ias
integer nr,ist
! allocatable local arrays
real(8), allocatable :: rfmt(:,:)
! external functions
real(8) rfmtinp
external rfmtinp
! allocate local arrays
allocate(rfmt(lmmaxvr,nrmtmax))
! calculate the kinetic energy for core states
engykncr=0.d0
do is=1,nspecies
  nr=nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! sum of core eigenvalues
    do ist=1,spnst(is)
      if (spcore(ist,is)) engykncr=engykncr+spocc(ist,is)*evalcr(ist,ias)
    end do
! core density
    if (spincore) then
      rfmt(1,1:nr)=(rhocr(1:nr,ias,1)+rhocr(1:nr,ias,2))/y00
    else
      rfmt(1,1:nr)=rhocr(1:nr,ias,1)/y00
    end if
    engykncr=engykncr-rfmtinp(1,0,nr,spr(:,is),lmmaxvr,rfmt,veffmt(:,:,ias))
  end do
end do
deallocate(rfmt)
return
end subroutine

