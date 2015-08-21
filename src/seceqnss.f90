
! Copyright (C) 2006 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine seceqnss(ik,apwalm,evalfv,evecfv,evecsv)
use modmain
use modldapu
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
real(8), intent(in) :: evalfv(nstfv,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ispn,jspn,is,ia,ias
integer ist,jst,i,j,k,l,lm,nm
integer nrc,ir,igk,ifg
integer lwork,info
real(8) cb,t1
real(8) ts0,ts1
complex(8) zq
! allocatable arrays
real(8), allocatable :: bir(:,:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: wfmt1(:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: wfmt3(:,:)
complex(8), allocatable :: wfmt4(:,:,:)
complex(8), allocatable :: wfir1(:,:)
complex(8), allocatable :: wfir2(:)
complex(8), allocatable :: zv(:,:)
complex(8), allocatable :: work(:)
! external functions
complex(8) zdotc,zfmtinp
external zdotc,zfmtinp
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(seceqnss): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
call timesec(ts0)
! coupling constant of the external field (g_e/4c)
cb=gfacte/(4.d0*solsc)
! zero the second-variational Hamiltonian (stored in the eigenvector array)
evecsv(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmaxvr,nrcmtmax,nspnfv))
allocate(wfmt3(lmmaxvr,nrcmtmax))
allocate(wfmt4(lmmaxvr,nrcmtmax,3))
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! de-phasing factor (FC, FB & LN)
    t1=-0.5d0*dot_product(vqcss(:),atposc(:,ia,is))
    zq=cmplx(cos(t1),sin(t1),8)
! compute the first-variational wavefunctions
    do ispn=1,nspnfv
      if (ispn.eq.2) zq=conjg(zq)
      do ist=1,nstfv
        call wavefmt(lradstp,lmaxvr,is,ia,ngk(ispn,ik),apwalm(:,:,:,:,ispn), &
         evecfv(:,ist,ispn),lmmaxvr,wfmt1(:,:,ist,ispn))
! de-phase if required
        if (ssdph) wfmt1(:,1:nrc,ist,ispn)=zq*wfmt1(:,1:nrc,ist,ispn)
      end do
    end do
    do jst=1,nstfv
      do ispn=1,nspnfv
! convert wavefunctions to spherical coordinates
        call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zbshtvr,lmmaxvr, &
         wfmt1(:,:,jst,ispn),lmmaxvr,zzero,wfmt2(:,:,ispn),lmmaxvr)
      end do
! apply effective magnetic field and convert to spherical harmonics
      wfmt3(:,1:nrc)=wfmt2(:,1:nrc,1)*beffmt(:,1:nrc,ias,3)
      call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr,wfmt3, &
       lmmaxvr,zzero,wfmt4(:,:,1),lmmaxvr)
      wfmt3(:,1:nrc)=-wfmt2(:,1:nrc,2)*beffmt(:,1:nrc,ias,3)
      call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr,wfmt3, &
       lmmaxvr,zzero,wfmt4(:,:,2),lmmaxvr)
      wfmt3(:,1:nrc)=wfmt2(:,1:nrc,2) &
       *cmplx(beffmt(:,1:nrc,ias,1),-beffmt(:,1:nrc,ias,2),8)
      call zgemm('N','N',lmmaxvr,nrc,lmmaxvr,zone,zfshtvr,lmmaxvr,wfmt3, &
       lmmaxvr,zzero,wfmt4(:,:,3),lmmaxvr)
! apply LDA+U potential if required
      if ((ldapu.ne.0).and.(llu(is).ge.0)) then
        l=llu(is)
        nm=2*l+1
        lm=idxlm(l,-l)
        do k=1,3
          if (k.eq.1) then
            ispn=1
            jspn=1
          else if (k.eq.2) then
            ispn=2
            jspn=2
          else
            ispn=1
            jspn=2
          end if
          call zgemm('N','N',nm,nrc,nm,zone,vmatlu(lm,lm,ispn,jspn,ias), &
           lmmaxlu,wfmt1(lm,1,jst,jspn),lmmaxvr,zone,wfmt4(lm,1,k),lmmaxvr)
        end do
      end if
! second-variational Hamiltonian matrix
      do ist=1,nstfv
        do k=1,3
          if (k.eq.1) then
            ispn=1
            i=ist
            j=jst
          else if (k.eq.2) then
            ispn=2
            i=ist+nstfv
            j=jst+nstfv
          else
            ispn=1
            i=ist
            j=jst+nstfv
          end if
          if (i.le.j) then
            evecsv(i,j)=evecsv(i,j)+zfmtinp(.true.,lmaxvr,nrc,rcmt(:,is), &
             lmmaxvr,wfmt1(:,:,ist,ispn),wfmt4(:,:,k))
          end if
        end do
      end do
    end do
! end loops over atoms and species
  end do
end do
deallocate(wfmt1,wfmt2,wfmt3,wfmt4)
!---------------------------!
!     interstitial part     !
!---------------------------!
allocate(bir(ngrtot,3))
allocate(wfir1(ngrtot,nspnfv))
allocate(wfir2(ngrtot))
allocate(zv(ngkmax,3))
do ir=1,ngrtot
  bir(ir,:)=(bxcir(ir,:)+cb*bfieldc(:))*cfunir(ir)
end do
do jst=1,nstfv
  do ispn=1,nspnfv
    wfir1(:,ispn)=0.d0
    do igk=1,ngk(ispn,ik)
      ifg=igfft(igkig(igk,ispn,ik))
      wfir1(ifg,ispn)=evecfv(igk,jst,ispn)
    end do
! Fourier transform wavefunction to real-space
    call zfftifc(3,ngrid,1,wfir1(:,ispn))
  end do
! multiply with magnetic field and transform to G-space
  do k=1,3
    if (k.eq.1) then
      ispn=1
      wfir2(:)=wfir1(:,1)*bir(:,3)
    else if (k.eq.2) then
      ispn=2
      wfir2(:)=-wfir1(:,2)*bir(:,3)
    else
      ispn=1
      wfir2(:)=wfir1(:,2)*cmplx(bir(:,1),-bir(:,2),8)
    end if
    call zfftifc(3,ngrid,-1,wfir2)
    do igk=1,ngk(ispn,ik)
      ifg=igfft(igkig(igk,ispn,ik))
      zv(igk,k)=wfir2(ifg)
    end do
  end do
  do ist=1,nstfv
    do k=1,3
      if (k.eq.1) then
        ispn=1
        i=ist
        j=jst
      else if (k.eq.2) then
        ispn=2
        i=ist+nstfv
        j=jst+nstfv
      else
        ispn=1
        i=ist
        j=jst+nstfv
      end if
      if (i.le.j) then
        evecsv(i,j)=evecsv(i,j) &
         +zdotc(ngk(ispn,ik),evecfv(:,ist,ispn),1,zv(:,k),1)
      end if
    end do
  end do
end do
deallocate(bir,wfir1,wfir2,zv)
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist,ispn)
  end do
end do
! diagonalise the second-variational Hamiltonian
allocate(rwork(3*nstsv))
lwork=2*nstsv
allocate(work(lwork))
call zheev('V','U',nstsv,evecsv,nstsv,evalsv(:,ik),work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(seceqnss): diagonalisation of the second-variational &
   &Hamiltonian failed")')
  write(*,'(" for k-point ",I8)') ik
  write(*,'(" ZHEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
deallocate(rwork,work)
call timesec(ts1)
!$OMP CRITICAL
timesv=timesv+ts1-ts0
!$OMP END CRITICAL
return
end subroutine

