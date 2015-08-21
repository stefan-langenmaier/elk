
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writewfpw
use modmain
implicit none
! local variables
integer ik,ist,igk,recl
real(8) chg,wfn0,wfn1,sum
! allocatable arrays
complex(8), allocatable :: wfpw(:,:,:)
complex(8), allocatable :: wfpwh(:,:,:,:,:)
! initialise global variables
call init0
call init1
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! determine the record length and open WFPW.OUT
allocate(wfpw(ngkmax,nspinor,nstsv))
inquire(iolength=recl) vkl(:,ik),ngkmax,nspinor,nstsv,wfpw
deallocate(wfpw)
open(50,file='WFPW.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 status='REPLACE',recl=recl)
! determine the record length and open WFPWH.OUT
allocate(wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
inquire(iolength=recl) vkl(:,ik),ngkmax,nspinor,nstsv,wfpwh
deallocate(wfpwh)
open(51,file='WFPWH.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
 status='REPLACE',recl=recl)
chg=0.d0
wfn0=1.d0
wfn1=0.d0
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfpw,wfpwh) &
!$OMP PRIVATE(ist,sum,igk)
!$OMP DO
do ik=1,nkpt
  allocate(wfpw(ngkmax,nspinor,nstsv))
  allocate(wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
!$OMP CRITICAL
  write(*,'("Info(writewfpw): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
  call genwfpw(vkl(:,ik),ngk(1,ik),igkig(:,1,ik),vgkl(:,:,1,ik),gkc(:,1,ik), &
   tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),wfpw,wfpwh)
!$OMP CRITICAL
  write(50,rec=ik) vkl(:,ik),ngkmax,nspinor,nstsv,wfpw
  write(51,rec=ik) vkl(:,ik),ngkmax,nspinor,nstsv,wfpwh
!$OMP END CRITICAL
! read in the occupancies
  call getoccsv(vkl(:,ik),occsv(:,ik))
! calculate the total charge from the low G+k wavefunctions
  do ist=1,nstsv
    sum=0.d0
    do igk=1,ngk(1,ik)
      sum=sum+dble(wfpw(igk,1,ist))**2+aimag(wfpw(igk,1,ist))**2
      if (spinpol) sum=sum+dble(wfpw(igk,2,ist))**2+aimag(wfpw(igk,2,ist))**2
    end do
    if (sum.lt.wfn0) wfn0=sum
    if (sum.gt.wfn1) wfn1=sum
!$OMP CRITICAL
    chg=chg+wkpt(ik)*occsv(ist,ik)*sum
!$OMP END CRITICAL
  end do
  deallocate(wfpw,wfpwh)
end do
!$OMP END DO
!$OMP END PARALLEL
close(50)
close(51)
write(*,*)
write(*,'("Info(writewfpw): low and high plane wave wavefunctions written to")')
write(*,'(" WFPW.OUT and WFPWH.OUT, respectively")')
write(*,*)
write(*,'(" Wavefunction norms")')
write(*,'("  minimum : ",G18.10)') wfn0
write(*,'("  maximum : ",G18.10)') wfn1
write(*,*)
write(*,'(" Total charge from low plane waves : ",G18.10)') chg
write(*,'(" Total valence charge              : ",G18.10)') chgval
write(*,*)
return
end subroutine

