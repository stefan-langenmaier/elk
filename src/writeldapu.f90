
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writeldapu
use modmain
use modldapu
implicit none
! local variables
integer is,ia,ias,ispn,jspn
integer l,k,m1,m2,lm1,lm2
! machine readable density matrix file
open(50,file='DMATLU'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  l=llu(is)
  if (l.ge.0) then
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,*)
      write(50,*)
      write(50,'(3I4," : species, atom, l")') is,ia,l
      do ispn=1,nspinor
        do jspn=1,nspinor
          write(50,*)
          write(50,'(2I4," : ispn, jspn; m1, m2, dmatlu below")') ispn,jspn
          do m1=-l,l
            lm1=idxlm(l,m1)
            do m2=-l,l
              lm2=idxlm(l,m2)
              write(50,'(2I6," ",2G18.10)') m1,m2,dmatlu(lm1,lm2,ispn,jspn,ias)
            end do
          end do
        end do
      end do
    end do
  end if
end do
close(50)
! machine readable potential matrix file
open(50,file='VMATLU'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  l=llu(is)
  if (l.ge.0) then
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,*)
      write(50,*)
      write(50,'(3I4," : species, atom, l")') is,ia,l
      do ispn=1,nspinor
        do jspn=1,nspinor
          write(50,*)
          write(50,'(2I4," : ispn, jspn; m1, m2, vmatlu below")') ispn,jspn
          do m1=-l,l
            lm1=idxlm(l,m1)
            do m2=-l,l
              lm2=idxlm(l,m2)
              write(50,'(2I6," ",2G18.10)') m1,m2,vmatlu(lm1,lm2,ispn,jspn,ias)
            end do
          end do
        end do
      end do
    end do
  end if
end do
close(50)
! machine readable alpha parameters
if ((ldapu.eq.3).and.(.not.readalu)) then
  open(50,file='ALPHALU'//trim(filext),action='WRITE',form='FORMATTED')
  do is=1,nspecies
    l=llu(is)
    if (l.ge.0) then
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        write(50,*)
        write(50,'(3I4," : species, atom, l")') is,ia,l
        write(50,'(G18.10," : alpha")') alphalu(ias)
      end do
    end if
  end do
  close(50)
end if
! Slater parameters
open(50,file='FLU'//trim(filext),action='WRITE',form='FORMATTED')
do is=1,nspecies
  l=llu(is)
  if (l.ge.0) then
    do ia=1,natoms(is)
      write(50,*)
      write(50,'(3I4," : species, atom, l")') is,ia,l
      do k=0,2*l,2
        write(50,'(G18.10," : F^(",I1,")")') flu(k,is),k
      end do
      write(50,'(G18.10," : U")') ujlu(1,is)
      write(50,'(G18.10," : J")') ujlu(2,is)
      if (inptypelu.ge.4) write(50,'(G18.10," : lambdalu")') lambdalu(is)
    end do
  end if
end do
close(50)
return
end subroutine

