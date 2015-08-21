
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepresk(ik,vclcv,vclvv,dvxmt,dvxir,dbxmt,dbxir)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: vclcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(in) :: vclvv(nstsv,nstsv,nkpt)
real(8), intent(inout) :: dvxmt(lmmaxvr,nrcmtmax,natmtot)
real(8), intent(inout) :: dvxir(ngtot)
real(8), intent(inout) :: dbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
real(8), intent(inout) :: dbxir(ngtot,ndmag)
! local variables
integer ist,jst,is,ia,ias
integer nrc,nrci,ic,m,idm
real(8) de
complex(8) z1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:),wfcr(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zmagmt(:,:,:,:),zmagir(:,:)
complex(8), allocatable :: zvfmt(:,:,:)
! external functions
complex(8) zfinp,zfmtinp
external zfinp,zfmtinp
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngtot,nspinor,nstsv))
allocate(wfcr(lmmaxvr,nrcmtmax,2))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngtot))
if (spinpol) then
  allocate(zmagmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
  allocate(zmagir(ngtot,ndmag))
  allocate(zvfmt(lmmaxvr,nrcmtmax,ndmag))
end if
! get the eigenvalues/vectors from file for input k-point
call getevalsv(vkl(:,ik),evalsv(:,ik))
call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! calculate the wavefunctions for all states
call genwfsv(.false.,.false.,nstsv,idx,ngk(1,ik),igkig(:,1,ik),apwalm,evecfv, &
 evecsv,wfmt,ngtot,wfir)
!-----------------------------------------------------------!
!     core-conduction overlap density and magnetisation     !
!-----------------------------------------------------------!
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    ic=0
    do ist=1,spnst(is)
      if (spcore(ist,is)) then
        do m=-spk(ist,is),spk(ist,is)-1
          ic=ic+1
! pass in m-1/2 to wavefcr
          call wavefcr(.false.,lradstp,is,ia,ist,m,nrcmtmax,wfcr)
          do jst=1,nstsv
            if (evalsv(jst,ik).gt.efermi) then
! calculate the complex overlap density in the muffin-tin
              zrhomt(:,1:nrc,ias)=conjg(wfcr(:,1:nrc,1))*wfmt(:,1:nrc,ias,1,jst)
              if (spinpol) then
                zrhomt(:,1:nrc,ias)=zrhomt(:,1:nrc,ias) &
                 +conjg(wfcr(:,1:nrc,2))*wfmt(:,1:nrc,ias,2,jst)
              end if
              z1=conjg(vclcv(ic,ias,jst,ik))
              z1=z1-zfmtinp(.false.,nrc,nrci,rcmt(:,is),zrhomt(:,:,ias), &
               zvxmt(:,:,ias))
! spin-polarised case
              if (spinpol) then
                call genzmagmt(is,wfcr(:,:,1),wfcr(:,:,2),wfmt(:,:,ias,1,jst), &
                 wfmt(:,:,ias,2,jst),1,zvfmt)
! integral of magnetisation dot exchange field
                do idm=1,ndmag
                  z1=z1-zfmtinp(.false.,nrc,nrci,rcmt(:,is),zvfmt(:,:,idm), &
                   zbxmt(:,:,ias,idm))
                end do
! end spin-polarised case
              end if
              de=evalcr(ist,ias)-evalsv(jst,ik)
              z1=z1*occmax*wkpt(ik)/(de+zi*swidth)
! residuals for exchange potential and field
!$OMP CRITICAL
              dvxmt(:,1:nrc,ias)=dvxmt(:,1:nrc,ias)+dble(z1*zrhomt(:,1:nrc,ias))
              do idm=1,ndmag
                dbxmt(:,1:nrc,ias,idm)=dbxmt(:,1:nrc,ias,idm) &
                 +dble(z1*zvfmt(:,1:nrc,idm))
              end do
!$OMP END CRITICAL
! end loop over jst
            end if
          end do
        end do
! end loop over ist
      end if
    end do
! end loops over atoms and species
  end do
end do
!--------------------------------------------------------------!
!     valence-conduction overlap density and magnetisation     !
!--------------------------------------------------------------!
do ist=1,nstsv
  if (evalsv(ist,ik).lt.efermi) then
    do jst=1,nstsv
      if (evalsv(jst,ik).gt.efermi) then
! calculate the overlap density
        call genzrho(.false.,.true.,wfmt(:,:,:,:,ist),wfir(:,:,ist), &
         wfmt(:,:,:,:,jst),wfir(:,:,jst),zrhomt,zrhoir)
        z1=conjg(vclvv(ist,jst,ik))
        z1=z1-zfinp(.false.,zrhomt,zrhoir,zvxmt,zvxir)
! spin-polarised case
        if (spinpol) then
          call genzmag(wfmt(:,:,:,:,ist),wfmt(:,:,:,:,jst),wfir(:,:,ist), &
           wfir(:,:,jst),zmagmt,zmagir)
! integral of magnetisation dot exchange field
          do idm=1,ndmag
            z1=z1-zfinp(.false.,zmagmt(:,:,:,idm),zmagir(:,idm), &
             zbxmt(:,:,:,idm),zbxir(:,idm))
          end do
        end if
        de=evalsv(ist,ik)-evalsv(jst,ik)
        z1=z1*occmax*wkpt(ik)/(de+zi*swidth)
! residuals for exchange potential and field
!$OMP CRITICAL
        do is=1,nspecies
          nrc=nrcmt(is)
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            dvxmt(:,1:nrc,ias)=dvxmt(:,1:nrc,ias)+dble(z1*zrhomt(:,1:nrc,ias))
            do idm=1,ndmag
              dbxmt(:,1:nrc,ias,idm)=dbxmt(:,1:nrc,ias,idm) &
               +dble(z1*zmagmt(:,1:nrc,ias,idm))
            end do
          end do
        end do
        dvxir(:)=dvxir(:)+dble(z1*zrhoir(:))
        do idm=1,ndmag
          dbxir(:,idm)=dbxir(:,idm)+dble(z1*zmagir(:,idm))
        end do
!$OMP END CRITICAL
! end loop over jst
      end if
    end do
! end loop over ist
  end if
end do
deallocate(apwalm,evecfv,evecsv)
deallocate(wfmt,wfir,wfcr,zrhomt,zrhoir)
if (spinpol) then
  deallocate(zmagmt,zmagir,zvfmt)
end if
return
end subroutine

