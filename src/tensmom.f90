
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: tensmom
! !INTERFACE:
subroutine tensmom(fnum)
! !USES:
use modmain
use modldapu
use modtest
! !DESCRIPTION:
!   Decompose the density matrix and the LDA+U Hartree-Fock energy in tensor
!   moments components, see {\it Phys. Rev. B} {\bf 80}, 035121 (2009).
!
! !REVISION HISTORY:
!   Created April 2008 (F.Cricchio and L.Nordstrom)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer, parameter :: ldim=2*lmaxlu+1
integer is,ia,ias
integer ispn,jspn,i1,i2
integer l,m1,m2,lm1,lm2
integer p,k,r,rmin,t
real(8) hfe,edir,exch,t1,t2
real(8) dmpoltot,dmpol,dmpol0
! automatic arrays
complex(8) tmom(-ldim:ldim)
! allocatable arrays
real(8),allocatable :: hfexch(:,:,:,:)
real(8),allocatable :: tmom2(:,:,:,:)
complex(8),allocatable :: dmcomb(:,:)
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(tensmom): the current implementation is only valid for &
   &spinpol=.true.")')
  write(*,*)
  stop
end if
if (iscl.le.1) then
  write(fnum,*)
  write(fnum,'("Tensor moment decomposition of density matrix and &
   &Hartree-Fock energy")')
  write(fnum,'(" (see Phys. Rev. B. 80, 035121 (2009))")')
  write(fnum,'("  W.W = modulus square of tensor moment")')
  write(fnum,'("  Edir = direct energy term")')
  write(fnum,'("  Ex = exchange energy term")')
  write(fnum,'("  Pol = polarization of density matrix")')
end if
if (iscl.ge.1) then
  write(fnum,*)
  write(fnum,'("+-------------------------+")')
  write(fnum,'("| Iteration number : ",I4," |")') iscl
  write(fnum,'("+-------------------------+")')
end if
! allocate arrays
allocate(dmcomb(ldim*nspinor,ldim*nspinor))
allocate(tmom2(0:2*lmaxlu,0:1,0:(2*lmaxlu+1),natmtot))
allocate(hfexch(0:2*lmaxlu,0:1,0:(2*lmaxlu+1),natmtot))
tmom2(:,:,:,:)=0.d0
hfexch(:,:,:,:)=0.d0
dmpol0=0.d0
! begin loop over atoms and species
do is=1,nspecies
  l=llu(is)
  if (l.lt.0) goto 10
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! write info to TENSMOM.OUT
    write(fnum,*)
    write(fnum,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
     trim(spsymb(is)),ia
    write(fnum,'(" l = ",I2)') l
! zero density matrix with combined indices
    dmcomb(:,:)=0.d0
    i1=0
    do ispn=1,nspinor
      do m1=-l,l
        lm1=idxlm(l,m1)
        i1=i1+1
        i2=0
        do jspn=1,nspinor
          do m2=-l,l
            lm2=idxlm(l,m2)
            i2=i2+1
! density matrix with combined indeces
            dmcomb(i1,i2)=dmatlu(lm1,lm2,ispn,jspn,ias)
          end do
        end do
      end do
    end do
! zero total Hartree-Fock energy
    hfe=0.d0
! zero total polarization
    dmpoltot=0.d0
    write(fnum,*)
    do k=0,2*l
      do p=0,1
        rmin=abs(k-p)
        do r=rmin,k+p
! decompose density matrix in tensor moment components
          call dm2tensmom(l,nspinor,ldim,k,p,r,dmcomb,tmom)
          t1=0.d0
! calculate modulus square of tensor moment
          do t=-r,r
            t1=t1+(-1)**t*dble(tmom(-t)*tmom(t))
          end do
! energy components
          call tensmomengy(is,l,k,p,r,edir,exch)
! save square of tensmom modulus and exch energy to be written in test file
          tmom2(k,p,r,ias)=t1
          hfexch(k,p,r,ias)=exch*t1
! polarization terms
          call dmplz(l,k,p,r,t1,dmpol)
! write to file square of tensmom modulus, direct and exch energy, polarizations
          write(fnum,'("  k = ",I1,", p = ",I1,", r = ",I1)') k,p,r
          if ((k+p+r).eq.0) then
! for k,p,r = 0 save reference polarization
            dmpol0=dmpol
! for k,p,r = 0 do not write out the polarization
            write(fnum,'("   W.W =",F14.8,", Edir =",F14.8,", Ex =",F14.8)') &
             t1,edir*t1,exch*t1
          else
! relative polarization
            t2=dmpol/dmpol0
            write(fnum,'("   W.W =",F14.8,", Edir =",F14.8,", Ex =",F14.8, &
             &", Pol =",F14.8)') t1,edir*t1,exch*t1,t2
! total relative polarization (skipping 000 component)
            dmpoltot=dmpoltot+t2
          end if
          hfe=hfe+(edir+exch)*t1
          do t=-r,r
! write out single components of tensor moments
            write(fnum,'("    t = ",I2," : ",2F14.8)') t,tmom(t)
          end do
          write(fnum,*)
        end do
      end do
    end do
    write(fnum,*)
    write(fnum,'("  Total Hartree-Fock energy (without DC correction) : ", &
     &F14.8)') hfe
    write(fnum,'("  Total polarization of density matrix : ",F14.8)') dmpoltot
    write(fnum,*)
! end loop over atoms and species
  end do
10 continue
end do
! write test files if required
if (test) then
  t1=sqrt(sum(tmom2(:,:,:,:)**2))
  call writetest(820,'RMS of tensor moments',tol=1.d-4,rv=t1)
  t1=sqrt(sum(hfexch(:,:,:,:)**2))
  call writetest(830,'RMS of LDA+U Hartree-Fock exchange energies',tol=1.d-4, &
   rv=t1)
end if
deallocate(dmcomb,tmom2,hfexch)
return
end subroutine
!EOC

