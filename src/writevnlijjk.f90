
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writevnlijjk
! !INTERFACE:
subroutine writevnlijjk
! !USES:
use modmain
use modmpi
! !DESCRIPTION:
!   Generates non-local Coulomb matrix elements of the type $(i-jj-k)$ and
!   outputs them to the file {\tt VNLIJJK.OUT}. Also writes the real diagonal
!   of this matrix, $(i-jj-i)$, to {\tt VNLIJJI.OUT}.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,ist,recl
! allocatable arrays
real(8), allocatable :: vnlijji(:,:,:)
complex(8), allocatable :: vnlijjk(:,:,:,:)
! determine record length for vnlijji and open file
allocate(vnlijji(nstsv,nstsv,nkpt))
inquire(iolength=recl) vnlijji
deallocate(vnlijji)
open(100,file='VNLIJJI.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
! determine record length for vnlijjk and open file
allocate(vnlijjk(nstsv,nstsv,nstsv,nkpt))
inquire(iolength=recl) vnlijjk
deallocate(vnlijjk)
open(101,file='VNLIJJK.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 recl=recl)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(vnlijji,vnlijjk,ist)
!$OMP DO
do ik=1,nkptnr
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  allocate(vnlijji(nstsv,nstsv,nkpt))
  allocate(vnlijjk(nstsv,nstsv,nstsv,nkpt))
!$OMP CRITICAL
  write(*,'("Info(writevnlijjk): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL
! calculate non-local matrix elements of the type (i-jj-k)
  call genvnlijjk(ik,vnlijjk)
! make a copy of the diagonal elements (i-jj-i)
  do ist=1,nstsv
    vnlijji(ist,:,:)=dble(vnlijjk(ist,ist,:,:))
  end do
!$OMP CRITICAL
  write(100,rec=ik) vnlijji
  write(101,rec=ik) vnlijjk
!$OMP END CRITICAL
  deallocate(vnlijji,vnlijjk)
end do
!$OMP END DO
!$OMP END PARALLEL
close(100)
close(101)
return
end subroutine
!EOC
