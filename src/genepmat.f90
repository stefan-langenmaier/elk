
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genepmat(iq,vpl,dveffmt,dveffir,epmat)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
real(8), intent(in) :: vpl(3)
complex(8), intent(in) :: dveffmt(lmmaxvr,nrcmtmax,natmtot,3*natmtot)
complex(8), intent(in) :: dveffir(ngrtot,3*natmtot)
complex(8), intent(out) :: epmat(nstsv,nstsv,3*natmtot)
! local variables
integer lmax,lmmax
integer is,ia,ias
integer ngp,ngpq,igp
integer nrc,irc,ifg
integer ist,jst,ispn
integer i,j,k,l,m,n
real(8) vpc(3),vpql(3),vpqc(3)
complex(8) zt1
! allocatable arrays
integer, allocatable :: igpig(:)
integer, allocatable :: igpqig(:)
real(8), allocatable :: vgpl(:,:)
real(8), allocatable :: vgpc(:,:)
real(8), allocatable :: gpc(:)
real(8), allocatable :: tpgpc(:,:)
real(8), allocatable :: vgpql(:,:)
real(8), allocatable :: vgpqc(:,:)
real(8), allocatable :: gpqc(:)
real(8), allocatable :: tpgpqc(:,:)
complex(8), allocatable :: sfacgp(:,:)
complex(8), allocatable :: sfacgpq(:,:)
complex(8), allocatable :: apwalm1(:,:,:,:)
complex(8), allocatable :: apwalm2(:,:,:,:)
complex(8), allocatable :: evecfv1(:,:)
complex(8), allocatable :: evecfv2(:,:)
complex(8), allocatable :: evecsv1(:,:)
complex(8), allocatable :: evecsv2(:,:)
complex(8), allocatable :: wfmt1(:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: wfmt3(:,:)
complex(8), allocatable :: zfmt1(:,:)
complex(8), allocatable :: zfmt2(:,:)
complex(8), allocatable :: zfir1(:)
complex(8), allocatable :: zfir2(:)
complex(8), allocatable :: zv(:)
complex(8), allocatable :: epm(:,:,:)
! external functions
complex(8) zfmtinp,zdotc
external zfmtinp,zdotc
n=3*natmtot
lmax=min(lmaxmat,lmaxvr)
lmmax=(lmax+1)**2
! allocate local arrays
allocate(igpig(ngkmax))
allocate(igpqig(ngkmax))
allocate(vgpl(3,ngkmax))
allocate(vgpc(3,ngkmax))
allocate(gpc(ngkmax))
allocate(tpgpc(2,ngkmax))
allocate(vgpql(3,ngkmax))
allocate(vgpqc(3,ngkmax))
allocate(gpqc(ngkmax))
allocate(tpgpqc(2,ngkmax))
allocate(sfacgp(ngkmax,natmtot))
allocate(sfacgpq(ngkmax,natmtot))
allocate(apwalm1(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(apwalm2(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv1(nmatmax,nstfv))
allocate(evecfv2(nmatmax,nstfv))
if (tevecsv) then
  allocate(evecsv1(nstsv,nstsv))
  allocate(evecsv2(nstsv,nstsv))
end if
allocate(wfmt1(lmmaxvr,nrcmtmax))
allocate(wfmt2(lmmaxvr,nrcmtmax,nstfv))
allocate(wfmt3(lmmaxvr,nrcmtmax))
allocate(zfmt1(lmmaxvr,nrcmtmax))
allocate(zfmt2(lmmaxvr,nrcmtmax))
allocate(zfir1(ngrtot))
allocate(zfir2(ngrtot))
allocate(zv(ngkmax))
allocate(epm(nstfv,nstfv,n))
! p-vector in Cartesian coordinates
call r3mv(bvec,vpl,vpc)
! generate the G+p vectors
call gengpvec(vpl,vpc,ngp,igpig,vgpl,vgpc,gpc,tpgpc)
! generate the structure factors
call gensfacgp(ngp,vgpc,ngkmax,sfacgp)
! find the matching coefficients for k-point p
call match(ngp,gpc,tpgpc,sfacgp,apwalm1)
! get the eigenvectors for k-point p
call getevecfv(vpl,vgpl,evecfv1)
! p+q-vector in lattice coordinates
vpql(:)=vpl(:)+vql(:,iq)
! p+q-vector in Cartesian coordinates
call r3mv(bvec,vpql,vpqc)
! generate the G+p+q-vectors
call gengpvec(vpql,vpqc,ngpq,igpqig,vgpql,vgpqc,gpqc,tpgpqc)
! generate the structure factors
call gensfacgp(ngpq,vgpqc,ngkmax,sfacgpq)
! find the matching coefficients for k-point p+q
call match(ngpq,gpqc,tpgpqc,sfacgpq,apwalm2)
! get the eigenvectors for k-point p+q
call getevecfv(vpql,vgpql,evecfv2)
! set the first-variational matrix element array to zero
epm(:,:,:)=0.d0
!------------------------------------!
!     muffin-tin matrix elements     !
!------------------------------------!
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist=1,nstfv
! calculate the wavefunction for k-point p+q
      call wavefmt(lradstp,lmaxvr,is,ia,ngpq,apwalm2,evecfv2(:,ist),lmmaxvr, &
       wfmt1)
! convert from spherical harmonics to spherical coordinates
      call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zbshtvr,lmmaxvr,wfmt1, &
       lmmaxvr,zzero,wfmt2(:,:,ist),lmmaxvr)
    end do
    do jst=1,nstfv
! calculate the wavefunction for k-point p
      call wavefmt(lradstp,lmaxvr,is,ia,ngp,apwalm1,evecfv1(:,jst),lmmaxvr, &
       wfmt1)
! convert from spherical harmonics to spherical coordinates
      call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zbshtvr,lmmaxvr,wfmt1, &
       lmmaxvr,zzero,wfmt3,lmmaxvr)
      do ist=1,nstfv
! multiply wavefunctions together in real-space
        do irc=1,nrc
          zfmt1(:,irc)=conjg(wfmt3(:,irc))*wfmt2(:,irc,ist)
        end do
! convert from spherical coordinates to spherical harmonics
        call zgemm('N','N',lmmax,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr,zfmt1, &
         lmmaxvr,zzero,zfmt2,lmmaxvr)
        do i=1,n
          epm(ist,jst,i)=epm(ist,jst,i)+zfmtinp(.true.,lmax,nrc,rcmt(:,is), &
           lmmaxvr,zfmt2,dveffmt(:,:,ias,i))
        end do
      end do
    end do
! end loops over atoms and species
  end do
end do
!--------------------------------------!
!     interstitial matrix elements     !
!--------------------------------------!
! compute interstitial wavefunctions for k-point p
do jst=1,nstfv
  zfir1(:)=0.d0
  do igp=1,ngp
    ifg=igfft(igpig(igp))
    zfir1(ifg)=evecfv1(igp,jst)
  end do
! Fourier transform wavefunction to real-space
  call zfftifc(3,ngrid,1,zfir1)
! loop over phonon branches
  do i=1,n
! multiply the wavefunction with the change in effective potential
    zfir2(:)=zfir1(:)*dveffir(:,i)
! Fourier transform to G-space
    call zfftifc(3,ngrid,-1,zfir2)
! store as wavefunction with G+p+q index
    do igp=1,ngpq
      ifg=igfft(igpqig(igp))
      zv(igp)=zfir2(ifg)
    end do
! add to the first-variational matrix elements
    do ist=1,nstfv
      epm(ist,jst,i)=epm(ist,jst,i)+zdotc(ngpq,evecfv2(:,ist),1,zv,1)
    end do
  end do
end do
!-------------------------------------------!
!     second-variational matrix elements    !
!-------------------------------------------!
if (tevecsv) then
! get the second-variational eigenvectors
  call getevecsv(vpl,evecsv1)
  call getevecsv(vpql,evecsv2)
  epmat(:,:,:)=0.d0
  do i=1,nstsv
    do j=1,nstsv
      k=0
      do ispn=1,nspinor
        do ist=1,nstfv
          k=k+1
          l=(ispn-1)*nstfv
          do jst=1,nstfv
            l=l+1
            zt1=conjg(evecsv2(k,i))*evecsv1(l,j)
            do m=1,n
              epmat(i,j,m)=epmat(i,j,m)+epm(ist,jst,m)*zt1
            end do
          end do
        end do
      end do
    end do
  end do
else
  epmat(:,:,:)=epm(:,:,:)
end if
deallocate(igpig,igpqig,vgpl,vgpc,gpc,tpgpc,vgpql,vgpqc,gpqc,tpgpqc)
deallocate(sfacgp,sfacgpq,apwalm1,apwalm2,evecfv1,evecfv2)
if (tevecsv) deallocate(evecsv1,evecsv2)
deallocate(wfmt1,wfmt2,wfmt3,zfmt1,zfmt2,zfir1,zfir2,zv,epm)
return
end subroutine

