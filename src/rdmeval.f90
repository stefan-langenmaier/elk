
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmeval
! !INTERFACE:
subroutine rdmeval
! !USES:
use modrdm
use modmain
! !DESCRIPTION:
!   RDMFT eigenvalues are determined by calculating the derivative of the
!   total energy w.r.t. the occupation number at half the maximum occupancy
!   ($n_{\rm max}/2$). Outputs the eigenvalues and occupation numbers to the
!   file {\tt RDM\_EIGVAL.OUT}
!
! !REVISION HISTORY:
!   Created 2009 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,ist
real(8) t1
! allocatable arrays
real(8), allocatable :: dedn(:,:)
allocate(dedn(nstsv,nkpt))
do ik=1,nkpt
  do ist=1,nstsv
    t1=occsv(ist,ik)
    occsv(ist,ik)=occmax/2.d0
    call rdmdedn(dedn)
    evalsv(ist,ik)=-dedn(ist,ik)
    occsv(ist,ik)=t1
  end do
  call putevalsv(ik,evalsv(:,ik))
end do
deallocate(dedn)
! write out the RDMFT eigenvalues
open(99,file='RDM_EIGVAL'//trim(filext),action='WRITE',form='FORMATTED')
write(99,'(I6," : nkpt")') nkpt
write(99,'(I6," : nstsv")') nstsv
do ik=1,nkpt
  write(99,*)
  write(99,'(I6,3G18.10," : k-point, vkl")') ik,vkl(:,ik)
  write(99,'(" (state, eigenvalue and occupancy below)")')
  do ist=1,nstsv
    write(99,'(I6,2G18.10)') ist,evalsv(ist,ik),occsv(ist,ik)
  end do
  write(99,*)
end do
close(99)
return
end subroutine
!EOC
