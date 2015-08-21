
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: seceqnfv
! !INTERFACE:
subroutine seceqnfv(nmatp,ngp,igpig,vpc,vgpc,apwalm,evalfv,evecfv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
!   ngp    : number of G+k-vectors for augmented plane waves (in,integer)
!   igpig  : index from G+k-vectors to G-vectors (in,integer(ngkmax))
!   vpc    : k-vector in Cartesian coordinates (in,real(3))
!   vgpc   : G+k-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
! !DESCRIPTION:
!   Solves the secular equation,
!   $$ (H-\epsilon O)b=0, $$
!   for the all the first-variational states of the input $k$-point.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nmatp
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vpc(3)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer is,ia,np
real(8) v(1)
real(8) ts0,ts1
! allocatable arrays
complex(8), allocatable :: h(:),o(:)
if (tpmat) then
  np=(nmatp*(nmatp+1))/2
else
  np=nmatp**2
end if
!-----------------------------------------------!
!     Hamiltonian and overlap matrix set up     !
!-----------------------------------------------!
call timesec(ts0)
allocate(h(np),o(np))
!$OMP PARALLEL SECTIONS DEFAULT(SHARED) PRIVATE(is,ia)
!$OMP SECTION
! Hamiltonian
h(:)=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    call hmlaa(.false.,is,ia,ngp,apwalm,v,h)
    call hmlalo(.false.,is,ia,ngp,apwalm,v,h)
    call hmllolo(.false.,is,ia,ngp,v,h)
  end do
end do
call hmlistl(.false.,ngp,igpig,vgpc,v,h)
!$OMP SECTION
! overlap
o(:)=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    call olpaa(.false.,is,ia,ngp,apwalm,v,o)
    call olpalo(.false.,is,ia,ngp,apwalm,v,o)
    call olplolo(.false.,is,ia,ngp,v,o)
  end do
end do
call olpistl(.false.,ngp,igpig,v,o)
!$OMP END PARALLEL SECTIONS
call timesec(ts1)
!$OMP CRITICAL
timemat=timemat+ts1-ts0
!$OMP END CRITICAL
!------------------------------------!
!     solve the secular equation     !
!------------------------------------!
if (tseqr) then
! system has inversion symmetry: use real symmetric matrix eigen solver
  call seceqnfvr(nmatp,ngp,vpc,h,o,evalfv,evecfv)
else
! no inversion symmetry: use complex Hermititan matrix eigen solver
  call seceqnfvz(nmatp,h,o,evalfv,evecfv)
end if
deallocate(h,o)
return
end subroutine
!EOC

