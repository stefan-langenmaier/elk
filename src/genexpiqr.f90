
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genexpiqr(vpl,zqrmt,emat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
complex(8), intent(in) :: zqrmt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(out) :: emat(nstsv,nstsv)
! local variables
integer is,ia,ias
integer ngp,ngpq,igp,ifg
integer ist,jst,ispn
integer nrc,i,j,k,l
real(8) vpc(3),vpql(3),vpqc(3),t1
complex(8) zsum
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
complex(8), allocatable :: zfir(:)
complex(8), allocatable :: zv(:)
complex(8), allocatable :: em(:,:)
! external functions
real(8) gaunt
complex(8) zfmtinp,zdotc
external gaunt,zfmtinp,zdotc
! check if q-vector is zero
t1=abs(vecql(1))+abs(vecql(2))+abs(vecql(3))
if (t1.lt.epslat) then
  emat(:,:)=0.d0
  do i=1,nstsv
    emat(i,i)=1.d0
  end do
  return
end if
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
allocate(zfir(ngrtot))
allocate(zv(ngkmax))
allocate(em(nstfv,nstfv))
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
vpql(:)=vpl(:)+vecql(:)
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
em(:,:)=0.d0
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
! multiply by exp(-iq.r) (conjugate because zfmtinp conjugates first function)
      wfmt1(:,1:nrc)=wfmt2(:,1:nrc,ist)*conjg(zqrmt(:,1:nrc,ias))
! convert from spherical coordinates to spherical harmonics
      call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr,wfmt1, &
       lmmaxvr,zzero,wfmt2(:,:,ist),lmmaxvr)
    end do
    do jst=1,nstfv
! calculate the wavefunction for k-point p
      call wavefmt(lradstp,lmaxvr,is,ia,ngp,apwalm1,evecfv1(:,jst),lmmaxvr, &
       wfmt1)
      do ist=1,nstfv
        em(ist,jst)=em(ist,jst)+zfmtinp(.true.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
         wfmt2(:,:,ist),wfmt1)
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
  zfir(:)=0.d0
  do igp=1,ngp
    ifg=igfft(igpig(igp))
    zfir(ifg)=evecfv1(igp,jst)
  end do
! Fourier transform wavefunction to real-space
  call zfftifc(3,ngrid,1,zfir)
! multiply with the characteristic function
  zfir(:)=zfir(:)*cfunir(:)
! Fourier transform back to G-space
  call zfftifc(3,ngrid,-1,zfir)
! store as wavefunction with G+p+q index
  do igp=1,ngpq
    ifg=igfft(igpqig(igp))
    zv(igp)=zfir(ifg)
  end do
! add to the first-variational matrix elements
  do ist=1,nstfv
    em(ist,jst)=em(ist,jst)+zdotc(ngpq,evecfv2(:,ist),1,zv,1)
  end do
end do
!-------------------------------------------!
!     second-variational matrix elements    !
!-------------------------------------------!
if (tevecsv) then
! get the second-variational eigenvectors
  call getevecsv(vpl,evecsv1)
  call getevecsv(vpql,evecsv2)
  do i=1,nstsv
    do j=1,nstsv
      zsum=0.d0
      k=0
      do ispn=1,nspinor
        do ist=1,nstfv
          k=k+1
          l=(ispn-1)*nstfv
          do jst=1,nstfv
            l=l+1
            zsum=zsum+em(ist,jst)*conjg(evecsv2(k,i))*evecsv1(l,j)
          end do
        end do
      end do
      emat(i,j)=zsum
    end do
  end do
else
  emat(:,:)=em(:,:)
end if
deallocate(igpig,igpqig,vgpl,vgpc,gpc,tpgpc)
deallocate(vgpql,vgpqc,gpqc,tpgpqc,sfacgp,sfacgpq)
deallocate(apwalm1,apwalm2,evecfv1,evecfv2)
if (tevecsv) deallocate(evecsv1,evecsv2)
deallocate(wfmt1,wfmt2,zfir,zv,em)
return
end subroutine

