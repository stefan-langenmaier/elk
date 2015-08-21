
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genwfpw(vpl,ngp,igpig,vgpl,gpc,tpgpc,sfacgp,wfpw,wfpwh)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpl(3,ngkmax)
real(8), intent(in) :: gpc(ngkmax)
real(8), intent(in) :: tpgpc(2,ngkmax)
complex(8), intent(in) :: sfacgp(ngkmax,natmtot)
complex(8), intent(out) :: wfpw(ngkmax,nspinor,nstsv)
complex(8), intent(out) :: wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
! local variables
integer ist,ispn,igp,ifg
integer is,ia,ias,l,m,lm
integer nrc,irc
real(8) x,t1,t2,t3
complex(8) zsum,zt1,zt2
! automatic arrays
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
real(8) gr(nrcmtmax),cf(4,nrcmtmax)
complex(8) ylm(lmmaxvr)
! allocatable arrays
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
complex(8), allocatable :: zfft(:)
allocate(jl(0:lmaxvr,nrcmtmax))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor,nstsv))
allocate(zfft(ngrtot))
! get the eigenvectors
call getevecfv(vpl,vgpl,evecfv)
call getevecsv(vpl,evecsv)
! find the matching coefficients
call match(ngp,gpc,tpgpc,sfacgp,apwalm)
! calculate the second-variational wavefunctions for all states
call genwfsv(.true.,.false.,ngp,igpig,evalsv,apwalm,evecfv,evecsv,wfmt,wfir)
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
t1=sqrt(omega)
do ist=1,nstsv
  do ispn=1,nspinor
! multiply wavefunction by characteristic function
    zfft(:)=wfir(:,ispn,ist)*cfunir(:)
! Fourier transform to G-space
    call zfftifc(3,ngrid,-1,zfft)
    do igp=1,ngp
      ifg=igfft(igpig(igp))
      wfpw(igp,ispn,ist)=t1*zfft(ifg)
    end do
  end do
end do
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
! initialise the high G+p wavefunction
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    wfpwh(:,1:nrcmt(is),ias,:,:)=wfmt(:,1:nrcmt(is),ias,:,:)
  end do
end do
t1=fourpi/sqrt(omega)
! loop over (G+p)-vectors
do igp=1,ngp
  call genylm(lmaxvr,tpgpc(:,igp),ylm)
! loop over species
  do is=1,nspecies
    nrc=nrcmt(is)
! generate spherical Bessel functions
    do irc=1,nrc
      x=gpc(igp)*rcmt(irc,is)
      call sbessel(lmaxvr,x,jl(:,irc))
    end do
! loop over atoms
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! loop over states
      do ist=1,nstsv
        do ispn=1,nspinor
          do irc=1,nrc
            zsum=0.d0
            lm=0
            do l=0,lmaxvr
              zt1=jl(l,irc)*conjg(zil(l))
              do m=-l,l
                lm=lm+1
                zsum=zsum+zt1*wfmt(lm,irc,ias,ispn,ist)*ylm(lm)
              end do
            end do
            zsum=zsum*rcmt(irc,is)**2
            fr1(irc)=dble(zsum)
            fr2(irc)=aimag(zsum)
          end do
          call fderiv(-1,nrc,rcmt(:,is),fr1,gr,cf)
          t2=gr(nrc)
          call fderiv(-1,nrc,rcmt(:,is),fr2,gr,cf)
          t3=gr(nrc)
          zt2=t1*cmplx(t2,t3,8)
! low G+p wavefunction
          wfpw(igp,ispn,ist)=wfpw(igp,ispn,ist)+zt2*conjg(sfacgp(igp,ias))
! high G+p wavefunction
          do irc=1,nrc
            lm=0
            do l=0,lmaxvr
              zt1=t1*jl(l,irc)*zil(l)*zt2
              do m=-l,l
                lm=lm+1
                wfpwh(lm,irc,ias,ispn,ist)=wfpwh(lm,irc,ias,ispn,ist) &
                 -zt1*conjg(ylm(lm))
              end do
            end do
          end do
        end do
! end loop over states
      end do
! end loop over atoms
    end do
! end loop over species
  end do
! end loop over (G+p)-vectors
end do
deallocate(jl,apwalm,evecfv,evecsv)
deallocate(wfmt,wfir,zfft)
return
end subroutine

