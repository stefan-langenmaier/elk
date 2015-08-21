
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readspecies
use modmain
implicit none
! local variables
integer is,ist,iostat
integer io,nlx,ilx,lx,ilo
e0min=0.d0
do is=1,nspecies
  open(50,file=trim(sppath)//trim(spfname(is)),action='READ',status='OLD', &
   form='FORMATTED',iostat=iostat)
  if (iostat.ne.0) then
    write(*,*)
    write(*,'("Error(readspecies): error opening species file ",A)') &
     trim(sppath)//trim(spfname(is))
    write(*,*)
    stop
  end if
  read(50,*) spsymb(is)
  read(50,*) spname(is)
  read(50,*) spzn(is)
  read(50,*) spmass(is)
  read(50,*) sprmin(is),rmt(is),sprmax(is),nrmt(is)
  if (sprmin(is).le.0.d0) then
    write(*,*)
    write(*,'("Error(readspecies): sprmin <= 0 : ",G18.10)') sprmin(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (rmt(is).le.sprmin(is)) then
    write(*,*)
    write(*,'("Error(readspecies): rmt <= sprmin : ",2G18.10)') rmt(is), &
     sprmin(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (sprmax(is).lt.rmt(is)) then
    write(*,*)
    write(*,'("Error(readspecies): sprmax < rmt : ",2G18.10)') sprmax(is), &
     rmt(is)
    write(*,*)
    stop
  end if
  if (nrmt(is).lt.20) then
    write(*,*)
    write(*,'("Error(readspecies): nrmt too small : ",I8)') nrmt(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
! multiply nrmt by the scale factor
  nrmt(is)=nrmt(is)*nrmtscf
! reduce the minimum radial mesh point by the same factor
  sprmin(is)=sprmin(is)/dble(nrmtscf)
  read(50,*) spnst(is)
  if ((spnst(is).le.0).or.(spnst(is).gt.maxspst)) then
    write(*,*)
    write(*,'("Error(readspecies): spnst out of range : ",I8)') spnst(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  do ist=1,spnst(is)
    read(50,*) spn(ist,is),spl(ist,is),spk(ist,is),spocc(ist,is),spcore(ist,is)
    if (spn(ist,is).lt.1) then
      write(*,*)
      write(*,'("Error(readspecies): spn < 1 : ",I8)') spn(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (spl(ist,is).lt.0) then
      write(*,*)
      write(*,'("Error(readspecies): spl < 0 : ",I8)') spl(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (spk(ist,is).lt.1) then
      write(*,*)
      write(*,'("Error(readspecies): spk < 1 : ",I8)') spk(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (spocc(ist,is).lt.0.d0) then
      write(*,*)
      write(*,'("Error(readspecies): spocc < 0 : ",G18.10)') spocc(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
  end do
  read(50,*) apword(0,is)
  if (apword(0,is).le.0) then
    write(*,*)
    write(*,'("Error(readspecies): apword <= 0 : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (apword(0,is).gt.maxapword) then
    write(*,*)
    write(*,'("Error(readspecies): apword too large : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxapword in modmain and recompile code")')
    write(*,*)
    stop
  end if
! set the APW orders for l>0
  apword(1:lmaxapw,is)=apword(0,is)
  do io=1,apword(0,is)
    read(50,*) apwe0(io,0,is),apwdm(io,0,is),apwve(io,0,is)
    if (apwdm(io,0,is).lt.0) then
      write(*,*)
      write(*,'("Error(readspecies): apwdm < 0 : ",I8)') apwdm(io,0,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and order ",I4)') io
      write(*,*)
      stop
    end if
! set the APW linearisation energies, derivative orders and variability for l>0
    apwe0(io,1:lmaxapw,is)=apwe0(io,0,is)
    apwdm(io,1:lmaxapw,is)=apwdm(io,0,is)
    apwve(io,1:lmaxapw,is)=apwve(io,0,is)
    e0min=min(e0min,apwe0(io,0,is))
  end do
  read(50,*) nlx
  if (nlx.lt.0) then
    write(*,*)
    write(*,'("Error(readspecies): nlx < 0 : ",I8)') nlx
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  do ilx=1,nlx
    read(50,*) lx,io
    if (lx.lt.0) then
      write(*,*)
      write(*,'("Error(readspecies): lx < 0 : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    if (lx.gt.lmaxapw) then
      write(*,*)
      write(*,'("Error(readspecies): lx > lmaxapw : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    apword(lx,is)=io
    if (apword(lx,is).le.0) then
      write(*,*)
      write(*,'("Error(readspecies): apword <= 0 : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    if (apword(lx,is).gt.maxapword) then
      write(*,*)
      write(*,'("Error(readspecies): apword too large : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,'("Adjust maxapword in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do io=1,apword(lx,is)
      read(50,*) apwe0(io,lx,is),apwdm(io,lx,is),apwve(io,lx,is)
      if (apwdm(io,lx,is).lt.0) then
        write(*,*)
        write(*,'("Error(readspecies): apwdm < 0 : ",I8)') apwdm(io,lx,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" exception number ",I4)') ilx
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      end if
      e0min=min(e0min,apwe0(io,lx,is))
    end do
  end do
  read(50,*) nlorb(is)
  if (nlorb(is).lt.0) then
    write(*,*)
    write(*,'("Error(readspecies): nlorb < 0 : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (nlorb(is).gt.maxlorb) then
    write(*,*)
    write(*,'("Error(readspecies): nlorb too large : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxlorb in modmain and recompile code")')
    write(*,*)
    stop
  end if
  do ilo=1,nlorb(is)
    read(50,*) lorbl(ilo,is),lorbord(ilo,is)
    if (lorbl(ilo,is).lt.0) then
      write(*,*)
      write(*,'("Error(readspecies): lorbl < 0 : ",I8)') lorbl(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbl(ilo,is).gt.lmaxmat) then
      write(*,*)
      write(*,'("Error(readspecies): lorbl > lmaxmat : ",2I8)') lorbl(ilo,is), &
       lmaxmat
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbord(ilo,is).lt.2) then
      write(*,*)
      write(*,'("Error(readspecies): lorbord < 2 : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbord(ilo,is).gt.maxlorbord) then
      write(*,*)
      write(*,'("Error(readspecies): lorbord too large : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,'("Adjust maxlorbord in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do io=1,lorbord(ilo,is)
      read(50,*) lorbe0(io,ilo,is),lorbdm(io,ilo,is),lorbve(io,ilo,is)
      if (lorbdm(io,ilo,is).lt.0) then
        write(*,*)
        write(*,'("Error(readspecies): lorbdm < 0 : ",I8)') lorbdm(io,ilo,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" local-orbital ",I4)') ilo
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      end if
      e0min=min(e0min,lorbe0(io,ilo,is))
    end do
  end do
  close(50)
end do
! add conduction state local-orbitals if required
call addlorbcnd
return
end subroutine

