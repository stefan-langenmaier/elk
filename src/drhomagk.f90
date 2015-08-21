
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine drhomagk(vpl,ngp,ngpq,igpig,igpqig,vgpl,gpqc,tpgpqc,sfacgpq,dwfpw)
use modmain
use modphonon
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(in) :: ngp(nspnfv)
integer, intent(in) :: ngpq(nspnfv)
integer, intent(in) :: igpig(ngkmax,nspnfv)
integer, intent(in) :: igpqig(ngkmax,nspnfv)
real(8), intent(in) :: vgpl(3,ngkmax,nspnfv)
real(8), intent(in) :: gpqc(ngkmax,nspnfv)
real(8), intent(in) :: tpgpqc(2,ngkmax,nspnfv)
complex(8), intent(in) :: sfacgpq(ngkmax,natmtot,nspnfv)
complex(8), intent(in) :: dwfpw(ngkmax,nspinor,nstsv)
! local variables
integer ispn,jspn,jspn0,jspn1
integer ik,ist,isym,is,ia,ias
integer l,m,lm,itp,nrc,irc
integer igp,igpq,ifg,ir
real(8) wo,x,t1
complex(8) zt1,zt2,zt3,zt4
! automatic arrays
complex(8) ylm(lmmaxvr)
! allocatable arrays
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: wfpw(:,:,:)
complex(8), allocatable :: wfir(:,:),dwfir(:,:)
complex(8), allocatable :: dwfmt(:,:,:,:,:)
complex(8), allocatable :: wfpwh(:,:,:,:,:)
complex(8), allocatable :: zfmt1(:,:,:),zfmt2(:,:,:)
! find the k-point number
call findkpt(vpl,isym,ik)
!-----------------------------!
!     low plane wave part     !
!-----------------------------!
allocate(wfpw(ngkmax,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor))
allocate(dwfir(ngrtot,nspinor))
! read in the low plane wave part of the wavefunctions
call getwfpw(vpl,vgpl,wfpw)
do ist=1,nstsv
  wo=wkpt(ik)*occsv(ist,ik)
  if (abs(wo).gt.epsocc) then
    t1=wo/omega
    do ispn=1,nspinor
      if (spinsprl) then
        jspn=ispn
      else
        jspn=1
      end if
      wfir(:,ispn)=0.d0
      do igp=1,ngp(jspn)
        ifg=igfft(igpig(igp,jspn))
        wfir(ifg,ispn)=wfpw(igp,ispn,ist)
      end do
      call zfftifc(3,ngrid,1,wfir(:,ispn))
      dwfir(:,ispn)=0.d0
      do igpq=1,ngpq(jspn)
        ifg=igfft(igpqig(igpq,jspn))
        dwfir(ifg,ispn)=dwfpw(igpq,ispn,ist)
      end do
      call zfftifc(3,ngrid,1,dwfir(:,ispn))
    end do
! add to density and magnetisation derivative
    if (spinpol) then
! spin-polarised
      if (ncmag) then
! non-collinear
        do ir=1,ngrtot
          zt1=conjg(wfir(ir,1))*dwfir(ir,1)
          zt2=conjg(wfir(ir,2))*dwfir(ir,2)
          zt3=conjg(wfir(ir,1))*dwfir(ir,2)
          zt4=conjg(wfir(ir,2))*dwfir(ir,1)
          drhoir(ir)=drhoir(ir)+t1*(zt1+zt2)
          dmagir(ir,1)=dmagir(ir,1)+t1*(zt3+zt4)
          dmagir(ir,2)=dmagir(ir,2)+t1*zi*(-zt3+zt4)
          dmagir(ir,3)=dmagir(ir,3)+t1*(zt1-zt2)
        end do
      else
! collinear
        do ir=1,ngrtot
          zt1=conjg(wfir(ir,1))*dwfir(ir,1)
          zt2=conjg(wfir(ir,2))*dwfir(ir,2)
          drhoir(ir)=drhoir(ir)+t1*(zt1+zt2)
          dmagir(ir,1)=dmagir(ir,1)+t1*(zt1-zt2)
        end do
      end if
    else
! spin-unpolarised
      drhoir(:)=drhoir(:)+t1*conjg(wfir(:,1))*dwfir(:,1)
    end if
  end if
end do
deallocate(wfpw,wfir,dwfir)
!------------------------------!
!     high plane wave part     !
!------------------------------!
allocate(jl(0:lmaxvr,nrcmtmax))
allocate(dwfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
dwfmt(:,:,:,:,:)=0.d0
t1=fourpi/sqrt(omega)
! loop over first-variational spin components
do ispn=1,nspnfv
  if (spinsprl) then
    jspn0=ispn; jspn1=ispn
  else
    jspn0=1; jspn1=nspinor
  end if
! loop over G+p+q vectors
  do igpq=1,ngpq(ispn)
! generate the spherical harmonics Y_lm(G+p+q)
    call genylm(lmaxvr,tpgpqc(:,igpq,ispn),ylm)
! loop over species
    do is=1,nspecies
      nrc=nrcmt(is)
! generate spherical Bessel functions
      do irc=1,nrc
        x=gpqc(igpq,ispn)*rcmt(irc,is)
        call sbessel(lmaxvr,x,jl(:,irc))
      end do
! loop over atoms
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        zt1=t1*sfacgpq(igpq,ias,ispn)
        do irc=1,nrc
          lm=0
          do l=0,lmaxvr
            zt2=jl(l,irc)*zt1*zil(l)
            do m=-l,l
              lm=lm+1
              zt3=zt2*conjg(ylm(lm))
              do ist=1,nstsv
                do jspn=jspn0,jspn1
                  dwfmt(lm,irc,ias,jspn,ist)=dwfmt(lm,irc,ias,jspn,ist) &
                   +dwfpw(igpq,jspn,ist)*zt3
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do
deallocate(jl)
allocate(wfpwh(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(zfmt1(lmmaxvr,nrcmtmax,nspinor))
allocate(zfmt2(lmmaxvr,nrcmtmax,nspinor))
! read in the high plane wave part of the wavefunctions
call getwfpwh(vpl,wfpwh)
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist=1,nstsv
      wo=wkpt(ik)*occsv(ist,ik)
      do ispn=1,nspinor
! convert wavefunctions to spherical coordinates
        call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zbshtvr,lmmaxvr, &
         wfpwh(:,:,ias,ispn,ist),lmmaxvr,zzero,zfmt1(:,:,ispn),lmmaxvr)
        call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zbshtvr,lmmaxvr, &
         dwfmt(:,:,ias,ispn,ist),lmmaxvr,zzero,zfmt2(:,:,ispn),lmmaxvr)
      end do
      if (spinpol) then
! spin-polarised
        if (ncmag) then
! non-collinear
          do irc=1,nrc
            do itp=1,lmmaxvr
              zt1=conjg(zfmt1(itp,irc,1))*zfmt2(itp,irc,1)
              zt2=conjg(zfmt1(itp,irc,2))*zfmt2(itp,irc,2)
              zt3=conjg(zfmt1(itp,irc,1))*zfmt2(itp,irc,2)
              zt4=conjg(zfmt1(itp,irc,2))*zfmt2(itp,irc,1)
              drhomt(itp,irc,ias)=drhomt(itp,irc,ias)+wo*(zt1+zt2)
              dmagmt(itp,irc,ias,1)=dmagmt(itp,irc,ias,1)+wo*(zt3+zt4)
              dmagmt(itp,irc,ias,2)=dmagmt(itp,irc,ias,2)+wo*zi*(-zt3+zt4)
              dmagmt(itp,irc,ias,3)=dmagmt(itp,irc,ias,3)+wo*(zt1-zt2)
            end do
          end do
        else
! collinear
          do irc=1,nrc
            do itp=1,lmmaxvr
              zt1=conjg(zfmt1(itp,irc,1))*zfmt2(itp,irc,1)
              zt2=conjg(zfmt1(itp,irc,2))*zfmt2(itp,irc,2)
              drhomt(itp,irc,ias)=drhomt(itp,irc,ias)+wo*(zt1+zt2)
              dmagmt(itp,irc,ias,1)=dmagmt(itp,irc,ias,1)+wo*(zt1-zt2)
            end do
          end do
        end if
      else
! spin-unpolarised
        do irc=1,nrc
          drhomt(:,irc,ias)=drhomt(:,irc,ias) &
           +wo*conjg(zfmt1(:,irc,1))*zfmt2(:,irc,1)
        end do
      end if
    end do
  end do
end do
deallocate(dwfmt,wfpwh,zfmt1,zfmt2)
return
end subroutine

