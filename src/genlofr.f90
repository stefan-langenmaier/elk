
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genlofr
! !INTERFACE:
subroutine genlofr
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the local-orbital radial functions. This is done by integrating
!   the scalar relativistic Schr\"{o}dinger equation (or its energy deriatives)
!   at the current linearisation energies using the spherical part of the
!   effective potential. For each local-orbital, a linear combination of
!   {\tt lorbord} radial functions is constructed such that its radial
!   derivatives up to order ${\tt lorbord}-1$ are zero at the muffin-tin radius.
!   This function is normalised and the radial Hamiltonian applied to it. The
!   results are stored in the global array {\tt lofr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Copied to equivalent atoms, February 2010 (A. Kozhevnikov and JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,jas
integer np,nr,ir,ilo,io,jo
integer j,l,nn,info
real(8) t1
! automatic arrays
logical done(natmmax)
real(8) vr(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
real(8) p0(nrmtmax,maxlorbord),p1(nrmtmax)
real(8) q0(nrmtmax,maxlorbord),q1(nrmtmax,maxlorbord)
real(8) p0s(nrmtmax),q0s(nrmtmax),q1s(nrmtmax)
real(8) hp0(nrmtmax)
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: xa(:),ya(:)
real(8), allocatable :: a(:,:),b(:),c(:)
! external functions
real(8) polynom
external polynom
! polynomial order
np=max(maxlorbord+1,4)
allocate(ipiv(np))
allocate(xa(np),ya(np),c(np))
allocate(a(np,np),b(np))
do is=1,nspecies
  nr=nrmt(is)
  done(:)=.false.
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
! use spherical part of potential
      vr(1:nr)=veffmt(1,1:nr,ias)*y00
      do ilo=1,nlorb(is)
        l=lorbl(ilo,is)
        do jo=1,lorbord(ilo,is)
! integrate the radial Schrodinger equation
          call rschroddme(solsc,lorbdm(jo,ilo,is),l,0,lorbe(jo,ilo,ias),nprad, &
           nr,spr(:,is),vr,nn,p0(:,jo),p1,q0(:,jo),q1(:,jo))
! normalise radial functions
          fr(1:nr)=p0(1:nr,jo)**2
          call fderiv(-1,nr,spr(:,is),fr,gr,cf)
          t1=1.d0/sqrt(abs(gr(nr)))
          p0(1:nr,jo)=t1*p0(1:nr,jo)
          q0(1:nr,jo)=t1*q0(1:nr,jo)
          q1(1:nr,jo)=t1*q1(1:nr,jo)
! set up the matrix of radial derivatives
          do j=1,np
            ir=nr-np+j
            xa(j)=spr(ir,is)
            ya(j)=p0(ir,jo)/spr(ir,is)
          end do
          do io=1,lorbord(ilo,is)
            a(io,jo)=polynom(io-1,np,xa,ya,c,rmt(is))
          end do
        end do
! set up the target vector
        b(:)=0.d0
        b(lorbord(ilo,is))=1.d0
        call dgesv(lorbord(ilo,is),1,a,np,ipiv,b,np,info)
        if (info.ne.0) then
          write(*,*)
          write(*,'("Error(genlofr): degenerate local-orbital radial &
           &functions")')
          write(*,'(" for species ",I4)') is
          write(*,'(" atom ",I4)') ia
          write(*,'(" and local-orbital ",I4)') ilo
          write(*,'(" ZGESV returned INFO = ",I8)') info
          write(*,*)
          stop
        end if
! generate linear superposition of radial functions
        p0s(:)=0.d0
        q0s(:)=0.d0
        q1s(:)=0.d0
        do io=1,lorbord(ilo,is)
          t1=b(io)
          p0s(1:nr)=p0s(1:nr)+t1*p0(1:nr,io)
          q0s(1:nr)=q0s(1:nr)+t1*q0(1:nr,io)
          q1s(1:nr)=q1s(1:nr)+t1*q1(1:nr,io)
        end do
! normalise radial functions
        fr(1:nr)=p0s(1:nr)**2
        call fderiv(-1,nr,spr(:,is),fr,gr,cf)
        t1=1.d0/sqrt(abs(gr(nr)))
        p0s(1:nr)=t1*p0s(1:nr)
        q0s(1:nr)=t1*q0s(1:nr)
        q1s(1:nr)=t1*q1s(1:nr)
! apply the scalar relativistic Hamiltonian
        call rschrodapp(solsc,l,nr,spr(:,is),vr,p0s,q0s,q1s,hp0)
! divide by r and store in global array
        do ir=1,nr
          t1=1.d0/spr(ir,is)
          lofr(ir,1,ilo,ias)=t1*p0s(ir)
          lofr(ir,2,ilo,ias)=t1*hp0(ir)
        end do
      end do
      done(ia)=.true.
! copy to equivalent atoms
      do ja=1,natoms(is)
        if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
          jas=idxas(ja,is)
          do ilo=1,nlorb(is)
            lofr(1:nr,1,ilo,jas)=lofr(1:nr,1,ilo,ias)
            lofr(1:nr,2,ilo,jas)=lofr(1:nr,2,ilo,ias)
          end do
          done(ja)=.true.
        end if
      end do
    end if
  end do
end do
deallocate(ipiv,xa,ya,a,b,c)
return
end subroutine
!EOC
