
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: wavefmt
! !INTERFACE:
subroutine wavefmt(lrstp,lmax,is,ia,ngp,apwalm,evecfv,ld,wfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp  : radial step length (in,integer)
!   lmax   : maximum angular momentum required (in,integer)
!   is     : species number (in,integer)
!   ia     : atom number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv : first-variational eigenvector (in,complex(nmatmax))
!   ld     : leading dimension (in,integer)
!   wfmt   : complex muffin-tin wavefunction passed in as real array
!            (out,real(2,ld,*))
! !DESCRIPTION:
!   Calculates the first-variational wavefunction in the muffin-tin in terms of
!   a spherical harmonic expansion. For atom $\alpha$ and a particular $k$-point
!   ${\bf p}$, the $r$-dependent $(l,m)$-coefficients of the wavefunction for
!   the $i$th state are given by
!   $$ \Phi^{i{\bf p}}_{\alpha lm}(r)=\sum_{\bf G}b^{i{\bf p}}_{\bf G}
!    \sum_{j=1}^{M^{\alpha}_l}A^{\alpha}_{jlm}({\bf G+p})u^{\alpha}_{jl}(r)
!    +\sum_{j=1}^{N^{\alpha}}b^{i{\bf p}}_{(\alpha,j,m)}v^{\alpha}_j(r)
!    \delta_{l,l_j}, $$
!   where $b^{i{\bf p}}$ is the $i$th eigenvector returned from routine
!   {\tt seceqn}; $A^{\alpha}_{jlm}({\bf G+p})$ is the matching coefficient;
!   $M^{\alpha}_l$ is the order of the APW; $u^{\alpha}_{jl}$ is the APW radial
!   function; $N^{\alpha}$ is the number of local-orbitals; $v^{\alpha}_j$ is
!   the $j$th local-orbital radial function; and $(\alpha,j,m)$ is a compound
!   index for the location of the local-orbital in the eigenvector. See routines
!   {\tt genapwfr}, {\tt genlofr}, {\tt match} and {\tt seceqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Fixed description, October 2004 (C. Brouder)
!   Removed argument ist, November 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp
integer, intent(in) :: lmax
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax)
integer, intent(in) :: ld
real(8), intent(out) :: wfmt(2,ld,*)
! local variables
integer ias,l,m,lm,ld2
integer ir,nr,io,ilo
! values smaller than eps are taken to be zero
real(8), parameter :: eps=1.d-14
complex(8) zt1
! external functions
complex(8) zdotu
external zdotu
if (lmax.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(wavefmt): lmax > lmaxapw : ",I8)') lmax
  write(*,*)
  stop
end if
ias=idxas(ia,is)
! zero the wavefunction
nr=0
do ir=1,nrmt(is),lrstp
  nr=nr+1
  wfmt(:,:,nr)=0.d0
end do
ld2=ld*2
! APW functions
do l=0,lmax
  do m=-l,l
    lm=idxlm(l,m)
    do io=1,apword(l,is)
      zt1=zdotu(ngp,evecfv,1,apwalm(:,io,lm,ias),1)
      if (abs(dble(zt1)).gt.1.d-14) then
        call daxpy(nr,dble(zt1),apwfr(:,1,io,l,ias),lrstp,wfmt(1,lm,1),ld2)
      end if
      if (abs(aimag(zt1)).gt.1.d-14) then
        call daxpy(nr,aimag(zt1),apwfr(:,1,io,l,ias),lrstp,wfmt(2,lm,1),ld2)
      end if
    end do
  end do
end do
! local-orbital functions
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  if (l.le.lmax) then
    do m=-l,l
      lm=idxlm(l,m)
      zt1=evecfv(ngp+idxlo(lm,ilo,ias))
      if (abs(dble(zt1)).gt.1.d-14) then
        call daxpy(nr,dble(zt1),lofr(:,1,ilo,ias),lrstp,wfmt(1,lm,1),ld2)
      end if
      if (abs(aimag(zt1)).gt.1.d-14) then
        call daxpy(nr,aimag(zt1),lofr(:,1,ilo,ias),lrstp,wfmt(2,lm,1),ld2)
      end if
    end do
  end if
end do
return
end subroutine
!EOC

