
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengy
use modmain
use modmpi
implicit none
! local variables
integer is,ia,nrc,m1,m2
integer ik,ist,jst
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: wfcr1(:,:,:)
complex(8), allocatable :: wfcr2(:,:,:)
complex(8), allocatable :: zrhomt(:,:)
complex(8), allocatable :: zvclmt(:,:)
complex(8), allocatable :: zfmt(:,:)
! external functions
complex(8) zfmtinp
external zfmtinp
allocate(wfcr1(lmmaxvr,nrcmtmax,2))
allocate(wfcr2(lmmaxvr,nrcmtmax,2))
allocate(zrhomt(lmmaxvr,nrcmtmax))
allocate(zvclmt(lmmaxvr,nrcmtmax))
allocate(zfmt(lmmaxvr,nrcmtmax))
! zero the exchange energy
engyx=0.d0
!--------------------------------------------------!
!     val-val-val and val-cr-val contributions     !
!--------------------------------------------------!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
!$OMP CRITICAL
  write(*,'("Info(exxengy): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
  call exxengyk(ik)
end do
!$OMP END DO
!$OMP END PARALLEL
! add energies from each process and redistribute
call mpi_allreduce(mpi_in_place,engyx,1,mpi_double_precision,mpi_sum, &
 mpi_comm_world,ierror)
!-----------------------------------!
!    core-core-core contribution    !
!-----------------------------------!
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    do jst=1,spnst(is)
      if (spcore(jst,is)) then
        do m2=-spk(jst,is),spk(jst,is)-1
          call wavefcr(lradstp,is,ia,jst,m2,nrcmtmax,wfcr2)
          do ist=1,spnst(is)
            if (spcore(ist,is)) then
              do m1=-spk(ist,is),spk(ist,is)-1
                call wavefcr(lradstp,is,ia,ist,m1,nrcmtmax,wfcr1)
! calculate the complex overlap density
                zfmt(:,1:nrc)=conjg(wfcr1(:,1:nrc,1))*wfcr2(:,1:nrc,1) &
                 +conjg(wfcr1(:,1:nrc,2))*wfcr2(:,1:nrc,2)
                call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr, &
                 zfmt,lmmaxvr,zzero,zrhomt,lmmaxvr)
! calculate the Coulomb potential
                call zpotclmt(lmaxvr,nrc,rcmt(:,is),lmmaxvr,zrhomt,zvclmt)
                zt1=zfmtinp(.true.,lmaxvr,nrc,rcmt(:,is),lmmaxvr,zrhomt,zvclmt)
                engyx=engyx-0.5d0*dble(zt1)
              end do
! end loop over ist
            end if
          end do
        end do
! end loop over jst
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(wfcr1,wfcr2,zrhomt,zvclmt,zfmt)
return
end subroutine

