
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpmat
! !INTERFACE:
subroutine genpmat(ngp,igpig,vgpc,apwalm,evecfv,evecsv,pmat)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer(nspnfv))
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax,nspnfv))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax,nspnfv))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pmat   : momentum matrix elements (out,complex(3,nstsv,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ P_{ij}=\int d^3r\,\Psi_{i{\bf k}}^*({\bf r})\left(-i\nabla
!    +\frac{1}{4c^2}\left[\vec{\sigma}\times\nabla V_s({\bf r})\right]\right)
!    \Psi_{j{\bf k}}({\bf r}), $$
!   where $V_s$ is the Kohn-Sham effective potential. The second term in the
!   brackets is only calculated if spin-orbit coupling is enabled.
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!   Fixed bug found by Juergen Spitaler, September 2006 (JKD)
!   Added spin-orbit correction, July 2010 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp(nspnfv)
integer, intent(in) :: igpig(ngkmax,nspnfv)
real(8), intent(in) :: vgpc(3,ngkmax,nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: pmat(3,nstsv,nstsv)
! local variables
integer ispn,jspn,ist,jst
integer is,ia,ias,nrc,ir,irc
integer igp,ifg,itp,i
complex(8) z11,z12,z21,z22,z31,z32,zt1,zt2
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
complex(8), allocatable :: gwfmt(:,:,:,:)
complex(8), allocatable :: gwfir(:,:)
complex(8), allocatable :: gvmt(:,:,:)
complex(8), allocatable :: zfmt1(:,:,:)
complex(8), allocatable :: zfmt2(:,:,:,:)
complex(8), allocatable :: zv(:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngkmax,nspinor,nstsv))
! calculate the second-variational wavefunctions for all states
call genwfsv(.true.,.true.,.false.,ngp,igpig,evalsv,apwalm,evecfv,evecsv,wfmt, &
 ngkmax,wfir)
! zero the momentum matrix elements array
pmat(:,:,:)=0.d0
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
allocate(gwfmt(lmmaxvr,nrcmtmax,3,nspinor))
if (spinorb) then
  allocate(gvmt(lmmaxvr,nrcmtmax,3))
  allocate(zfmt1(lmmaxvr,nrcmtmax,nspinor))
  allocate(zfmt2(lmmaxvr,nrcmtmax,3,nspinor))
end if
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute gradient of potential for spin-orbit correction if required
    if (spinorb) then
      irc=0
      do ir=1,nrmt(is),lradstp
        irc=irc+1
        call rtozflm(lmaxvr,veffmt(:,ir,ias),zfmt1(:,irc,1))
      end do
      call gradzfmt(lmaxvr,nrc,rcmt(:,is),lmmaxvr,nrcmtmax,zfmt1,gvmt)
! convert to spherical coordinates
      do i=1,3
        zfmt1(:,1:nrc,1)=gvmt(:,1:nrc,i)
        call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zbshtvr,lmmaxvr,zfmt1, &
         lmmaxvr,zzero,gvmt(:,:,i),lmmaxvr)
      end do
    end if
    do jst=1,nstsv
      do ispn=1,nspinor
! compute the gradient of the wavefunction
        call gradzfmt(lmaxvr,nrc,rcmt(:,is),lmmaxvr,nrcmtmax, &
         wfmt(:,:,ias,ispn,jst),gwfmt(:,:,:,ispn))
      end do
! add spin-orbit correction if required
      if (spinorb) then
        do ispn=1,nspinor
! convert wavefunction to spherical coordinates
          call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zbshtvr,lmmaxvr, &
           wfmt(:,:,ias,ispn,jst),lmmaxvr,zzero,zfmt1(:,:,ispn),lmmaxvr)
        end do
! compute sigma x (grad V(r)) psi(r)
        do irc=1,nrc
          do itp=1,lmmaxvr
            zt1=zfmt1(itp,irc,1)
            zt2=zfmt1(itp,irc,2)
            z11=gvmt(itp,irc,1)*zt1
            z12=gvmt(itp,irc,1)*zt2
            z21=gvmt(itp,irc,2)*zt1
            z22=gvmt(itp,irc,2)*zt2
            z31=gvmt(itp,irc,3)*zt1
            z32=gvmt(itp,irc,3)*zt2
            zfmt2(itp,irc,1,1)=-z21-cmplx(-aimag(z32),dble(z32),8)
            zfmt2(itp,irc,1,2)=z22+cmplx(-aimag(z31),dble(z31),8)
            zfmt2(itp,irc,2,1)=z11-z32
            zfmt2(itp,irc,2,2)=-z12-z31
            zfmt2(itp,irc,3,1)=cmplx(-aimag(z12),dble(z12),8)+z22
            zfmt2(itp,irc,3,2)=-cmplx(-aimag(z11),dble(z11),8)+z21
          end do
        end do
! convert to spherical harmonics and add to wavefunction gradient
        zt1=1.d0/(4.d0*solsc**2)
        do ispn=1,nspinor
          do i=1,3
            call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zt1,zfshtvr,lmmaxvr, &
             zfmt2(:,:,i,ispn),lmmaxvr,zone,gwfmt(:,:,i,ispn),lmmaxvr)
          end do
        end do
      end if
! find the overlaps
      do ispn=1,nspinor
        do ist=1,jst
          do i=1,3
            pmat(i,ist,jst)=pmat(i,ist,jst)+zfmtinp(.true.,lmaxvr,nrc, &
             rcmt(:,is),lmmaxvr,wfmt(:,:,ias,ispn,ist),gwfmt(:,:,i,ispn))
          end do
        end do
      end do
    end do
  end do
end do
deallocate(gwfmt)
if (spinorb) deallocate(gvmt,zfmt1,zfmt2)
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
allocate(gwfir(ngrtot,3),zv(ngkmax))
do jst=1,nstsv
  do ispn=1,nspinor
    if (spinsprl) then
      jspn=ispn
    else
      jspn=1
    end if
! compute the gradient
    gwfir(:,:)=0.d0
    do igp=1,ngp(jspn)
      ifg=igfft(igpig(igp,jspn))
      zt1=wfir(igp,ispn,jst)
      gwfir(ifg,:)=vgpc(:,igp,jspn)*cmplx(-aimag(zt1),dble(zt1),8)
    end do
    do i=1,3
! Fourier transform to real-space
      call zfftifc(3,ngrid,1,gwfir(:,i))
! multiply by characteristic function
      gwfir(:,i)=gwfir(:,i)*cfunir(:)
! Fourier transform back to G-space
      call zfftifc(3,ngrid,-1,gwfir(:,i))
    end do
! find the overlaps
    do i=1,3
      do igp=1,ngp(jspn)
        ifg=igfft(igpig(igp,jspn))
        zv(igp)=gwfir(ifg,i)
      end do
      do ist=1,jst
        pmat(i,ist,jst)=pmat(i,ist,jst)+zdotc(ngp(jspn),wfir(:,ispn,ist),1,zv,1)
      end do
    end do
  end do
end do
deallocate(gwfir,zv)
! multiply by -i and set lower triangular part
do ist=1,nstsv
  do jst=ist,nstsv
    do i=1,3
      zt1=pmat(i,ist,jst)
      zt1=cmplx(aimag(zt1),-dble(zt1),8)
      pmat(i,ist,jst)=zt1
      pmat(i,jst,ist)=conjg(zt1)
    end do
  end do
end do
deallocate(wfmt,wfir)
return
end subroutine
!EOC

