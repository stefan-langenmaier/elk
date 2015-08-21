
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genwfpw(twfpwh,vpl,ngp,igpig,vgpl,gpc,tpgpc,sfacgp,wfpw,wfpwh)
use modmain
implicit none
! arguments
logical, intent(in) :: twfpwh
real(8), intent(in) :: vpl(3)
integer, intent(in) :: ngp(nspnfv)
integer, intent(in) :: igpig(ngkmax,nspnfv)
real(8), intent(in) :: vgpl(3,ngkmax,nspnfv)
real(8), intent(in) :: gpc(ngkmax,nspnfv)
real(8), intent(in) :: tpgpc(2,ngkmax,nspnfv)
complex(8), intent(in) :: sfacgp(ngkmax,natmtot,nspnfv)
complex(8), intent(out) :: wfpw(ngkmax,nspinor,nstsv)
complex(8), intent(out) :: wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
! local variables
integer ispn0,ispn1,ispn,jspn
integer ist,igp,ifg,nrc,irc
integer is,ia,ias,l,m,lm
real(8) x,t0,t1,t2
complex(8) zsum1,zsum2,zt1,zt2,zt3
! automatic arrays
real(8) fr1(nrcmtmax),fr2(nrcmtmax),gr(nrcmtmax)
complex(8) ylm(lmmaxvr)
! allocatable arrays
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: zfft(:)
allocate(jl(0:lmaxvr,nrcmtmax))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor,nstsv))
allocate(zfft(ngrtot))
! get the eigenvectors from file
call getevecfv(vpl,vgpl,evecfv)
call getevecsv(vpl,evecsv)
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngp(ispn),gpc(:,ispn),tpgpc(:,:,ispn),sfacgp(:,:,ispn), &
   apwalm(:,:,:,:,ispn))
end do
! calculate the second-variational wavefunctions for all states
call genwfsv(.true.,.false.,.false.,ngp,igpig,evalsv,apwalm,evecfv,evecsv, &
 wfmt,ngrtot,wfir)
!-----------------------------------!
!     interstitial contribution     !
!-----------------------------------!
t0=sqrt(omega)
do ist=1,nstsv
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! multiply wavefunction by characteristic function
    zfft(:)=wfir(:,ispn,ist)*cfunir(:)
! Fourier transform to G-space
    call zfftifc(3,ngrid,-1,zfft)
    do igp=1,ngp(jspn)
      ifg=igfft(igpig(igp,jspn))
      wfpw(igp,ispn,ist)=t0*zfft(ifg)
    end do
  end do
end do
!---------------------------------!
!     muffin-tin contribution     !
!---------------------------------!
! initialise the high G+p wavefunction if required
if (twfpwh) then
  do ias=1,natmtot
    is=idxis(ias)
    wfpwh(:,1:nrcmt(is),ias,:,:)=wfmt(:,1:nrcmt(is),ias,:,:)
  end do
end if
t0=fourpi/sqrt(omega)
! loop over first-variational spin components
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
! loop over G+p-vectors
  do igp=1,ngp(jspn)
! generate the spherical harmonics Y_lm(G+p)
    call genylm(lmaxvr,tpgpc(:,igp,jspn),ylm)
! loop over species
    do is=1,nspecies
      nrc=nrcmt(is)
! generate spherical Bessel functions
      do irc=1,nrc
        x=gpc(igp,jspn)*rcmt(irc,is)
        call sbessel(lmaxvr,x,jl(:,irc))
      end do
! loop over atoms
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        zt1=conjg(sfacgp(igp,ias,jspn))
! loop over states
        do ist=1,nstsv
          do ispn=ispn0,ispn1
            do irc=1,nrc
              zsum1=0.d0
              lm=0
              do l=0,lmaxvr
                zsum2=0.d0
                do m=-l,l
                  lm=lm+1
                  zsum2=zsum2+wfmt(lm,irc,ias,ispn,ist)*ylm(lm)
                end do
                zsum1=zsum1+jl(l,irc)*conjg(zil(l))*zsum2
              end do
              zsum1=zsum1*rcmt(irc,is)**2
              fr1(irc)=dble(zsum1)
              fr2(irc)=aimag(zsum1)
            end do
            call fderiv(-1,nrc,rcmt(:,is),fr1,gr)
            t1=gr(nrc)
            call fderiv(-1,nrc,rcmt(:,is),fr2,gr)
            t2=gr(nrc)
            zt2=t0*cmplx(t1,t2,8)
! low G+p wavefunction
            wfpw(igp,ispn,ist)=wfpw(igp,ispn,ist)+zt1*zt2
! high G+p wavefunction if required
            if (twfpwh) then
              do irc=1,nrc
                lm=0
                do l=0,lmaxvr
                  zt3=t0*jl(l,irc)*zt2*zil(l)
                  do m=-l,l
                    lm=lm+1
                    wfpwh(lm,irc,ias,ispn,ist)=wfpwh(lm,irc,ias,ispn,ist) &
                     -zt3*conjg(ylm(lm))
                  end do
                end do
              end do
            end if
          end do
! end loop over states
        end do
! end loop over atoms
      end do
! end loop over species
    end do
! end loop over G+p-vectors
  end do
! end loop over first-variational spin components
end do
deallocate(jl,apwalm,evecfv,evecsv)
deallocate(wfmt,wfir,zfft)
return
end subroutine

