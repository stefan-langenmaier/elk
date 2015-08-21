
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: wavefmt
! !INTERFACE:
subroutine wavefmt(lrstp,lmax,ias,ngp,apwalm,evecfv,ld,wfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp  : radial step length (in,integer)
!   lmax   : maximum angular momentum required (in,integer)
!   ias    : joint atom and species number (in,integer)
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
!   {\tt eveqn}; $A^{\alpha}_{jlm}({\bf G+p})$ is the matching coefficient;
!   $M^{\alpha}_l$ is the order of the APW; $u^{\alpha}_{jl}$ is the APW radial
!   function; $N^{\alpha}$ is the number of local-orbitals; $v^{\alpha}_j$ is
!   the $j$th local-orbital radial function; and $(\alpha,j,m)$ is a compound
!   index for the location of the local-orbital in the eigenvector. See routines
!   {\tt genapwfr}, {\tt genlofr}, {\tt match} and {\tt eveqn}.
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
integer, intent(in) :: ias
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax)
integer, intent(in) :: ld
real(8), intent(out) :: wfmt(2,ld,*)
! local variables
integer is,nrc,ir,ld2
integer l,m,lm,io,ilo
complex(8) z1
! external functions
complex(8) zdotu
external zdotu
if (lmax.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(wavefmt): lmax > lmaxapw : ",I8)') lmax
  write(*,*)
  stop
end if
ld2=ld*2
! species number
is=idxis(ias)
! zero the wavefunction
nrc=0
do ir=1,nrmt(is),lrstp
  nrc=nrc+1
  wfmt(:,:,nrc)=0.d0
end do
! APW functions
lm=0
do l=0,lmax
  do m=-l,l
    lm=lm+1
    do io=1,apword(l,is)
      z1=zdotu(ngp,evecfv,1,apwalm(:,io,lm,ias),1)
      if (abs(dble(z1)).gt.1.d-14) then
        call daxpy(nrc,dble(z1),apwfr(:,1,io,l,ias),lrstp,wfmt(1,lm,1),ld2)
      end if
      if (abs(aimag(z1)).gt.1.d-14) then
        call daxpy(nrc,aimag(z1),apwfr(:,1,io,l,ias),lrstp,wfmt(2,lm,1),ld2)
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
      z1=evecfv(ngp+idxlo(lm,ilo,ias))
      if (abs(dble(z1)).gt.1.d-14) then
        call daxpy(nrc,dble(z1),lofr(:,1,ilo,ias),lrstp,wfmt(1,lm,1),ld2)
      end if
      if (abs(aimag(z1)).gt.1.d-14) then
        call daxpy(nrc,aimag(z1),lofr(:,1,ilo,ias),lrstp,wfmt(2,lm,1),ld2)
      end if
    end do
  end if
end do
return
end subroutine
!EOC

