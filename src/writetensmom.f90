
! Copyright (C) 2009 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetensmom
use modmain
use modldapu
implicit none
if (ldapu.eq.0) then
  write(*,*)
  write(*,'("Error(writetensmom): ldapu = 0")')
  write(*,*)
  stop
end if
! initialize universal variables
call init0
! read density matrix from file DMATLU.OUT
call readdmatlu
! read Slater integrals from FLU.OUT
call readflu
! open TENSMOM.OUT
open(50,file='TENSMOM.OUT',action='WRITE',form='FORMATTED')
! write the tensor
call tensmom(50)
close(50)
write(*,*)
write(*,'("Info(writetensmom): Tensor moment decomposition of density matrix")')
write(*,'(" and Hartree-Fock energy written to TENSMOM.OUT")')
write(*,*)
return
end subroutine

