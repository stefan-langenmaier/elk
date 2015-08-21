
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readdmatlu
use modmain
use modldapu
implicit none
! local variables
integer is,ia,ias,ispn,jspn
integer is_,ia_,ispn_,jspn_
integer l,m1,m2,lm1,lm2
integer l_,m1_,m2_
real(8) a,b
! read density matrix from DMATLU.OUT
open(50,file='DMATLU'//trim(filext),action='READ',form='FORMATTED')
do is=1,nspecies
  l=llu(is)
  if (l.ge.0) then
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      read(50,*)
      read(50,*)
      read(50,*) is_,ia_,l_
      if ((is.ne.is_).or.(ia.ne.ia_).or.(l.ne.l_)) then
        write(*,*)
        write(*,'("Error(readdmatlu): differing is, ia or l")')
        write(*,'(" current    : ",3I4)') is,ia,l
        write(*,'(" DMATLU.OUT : ",3I4)') is_,ia_,l_
        write(*,*)
        stop
      end if
      do ispn=1,nspinor
        do jspn=1,nspinor
          read(50,*)
          read(50,*) ispn_,jspn_
          if ((ispn.ne.ispn_).or.(jspn.ne.jspn_)) then
            write(*,*)
            write(*,'("Error(readdmatlu): differing ispn or jspn")')
            write(*,'(" current    : ",2I4)') ispn,jspn
            write(*,'(" DMATLU.OUT : ",2I4)') ispn_,jspn_
            write(*,*)
            stop
          end if
          do m1=-l,l
            lm1=idxlm(l,m1)
            do m2=-l,l
              lm2=idxlm(l,m2)
              read(50,*) m1_,m2_,a,b
              dmatlu(lm1,lm2,ispn,jspn,ias)=cmplx(a,b,8)
              if ((m1.ne.m1_).or.(m2.ne.m2_)) then
                write(*,*)
                write(*,'("Error(readdmatlu): differing m1 or m2")')
                write(*,'(" current    : ",2I6)') m1,m2
                write(*,'(" DMATLU.OUT : ",2I6)') m1_,m2_
                write(*,*)
                stop
              end if
            end do
          end do
        end do
      end do
    end do
  end if
end do
close(50)
return
end subroutine

