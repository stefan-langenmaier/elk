
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl,
! F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnsv(ngp,igpig,vgpc,apwalm,evalfv,evecfv,evalsvp,evecsv)
use modmain
use moddftu
implicit none
! arguments
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(in) :: evalfv(nstfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
real(8), intent(out) :: evalsvp(nstsv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer nsc,nsd,ld,ist,jst
integer ispn,jspn,is,ias
integer nrc,nrci,iro,irc
integer l,lm,lmi,nm,igp,ifg
integer i,j,k,lwork,info
real(8) ca,t1
real(8) ts0,ts1
complex(8) z1
! automatic arrays
complex(8) zlflm(lmmaxvr,3)
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: wfmt1(:,:,:),wfmt2(:,:),wfmt3(:,:),wfmt4(:,:,:)
complex(8), allocatable :: gwfmt(:,:,:),wfir1(:),wfir2(:),z(:,:),work(:)
! external functions
complex(8) zdotc,zfmtinp
external zdotc,zfmtinp
! no calculation of second-variational eigenvectors
if (.not.tevecsv) then
  do i=1,nstsv
    evalsvp(i)=evalfv(i)
  end do
  evecsv(:,:)=0.d0
  do i=1,nstsv
    evecsv(i,i)=1.d0
  end do
  return
end if
call timesec(ts0)
! coupling constant of the external A-field (1/c)
ca=1.d0/solsc
! number of spin combinations after application of Hamiltonian
if (spinpol) then
  if (ncmag) then
    nsc=3
  else
    nsc=2
  end if
  nsd=2
else
  nsc=1
  nsd=1
end if
ld=lmmaxdm*nspinor
! zero the second-variational Hamiltonian (stored in the eigenvector array)
evecsv(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
lmi=lmmaxinr
! begin loop over atoms
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,wfmt2,wfmt3,wfmt4,gwfmt) &
!$OMP PRIVATE(is,nrc,nrci,iro,ist,jst,irc,zlflm) &
!$OMP PRIVATE(t1,l,nm,lm,i,j,k,ispn,jspn)
!$OMP DO
do ias=1,natmtot
  allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv),wfmt2(lmmaxvr,nrcmtmax))
  allocate(wfmt3(lmmaxvr,nrcmtmax),wfmt4(lmmaxvr,nrcmtmax,nsc))
  if (afieldpol) allocate(gwfmt(lmmaxvr,nrcmtmax,3))
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmtinr(is)
  iro=nrci+1
! compute the first-variational wavefunctions
  do ist=1,nstfv
    call wavefmt(lradstp,lmaxvr,ias,ngp,apwalm,evecfv(:,ist),lmmaxvr, &
     wfmt1(:,:,ist))
  end do
! begin loop over states
  do jst=1,nstfv
    if (spinpol) then
! convert wavefunction to spherical coordinates
      call zbsht(nrc,nrci,wfmt1(:,:,jst),wfmt2)
! apply Kohn-Sham effective magnetic field
      wfmt3(1:lmi,1:nrci)=bsmt(1:lmi,1:nrci,ias,ndmag)*wfmt2(1:lmi,1:nrci)
      wfmt3(:,iro:nrc)=bsmt(:,iro:nrc,ias,ndmag)*wfmt2(:,iro:nrc)
! convert to spherical harmonics and store in wfmt4
      call zfsht(nrc,nrci,wfmt3,wfmt4(:,:,1))
      wfmt4(1:lmi,1:nrci,2)=-wfmt4(1:lmi,1:nrci,1)
      wfmt4(:,iro:nrc,2)=-wfmt4(:,iro:nrc,1)
! non-collinear magnetic field
      if (ncmag) then
        wfmt3(1:lmi,1:nrci)=cmplx(bsmt(1:lmi,1:nrci,ias,1), &
         -bsmt(1:lmi,1:nrci,ias,2),8)*wfmt2(1:lmi,1:nrci)
        wfmt3(:,iro:nrc)=cmplx(bsmt(:,iro:nrc,ias,1), &
         -bsmt(:,iro:nrc,ias,2),8)*wfmt2(:,iro:nrc)
        call zfsht(nrc,nrci,wfmt3,wfmt4(:,:,3))
      end if
! apply spin-orbit coupling if required
      if (spinorb) then
! inner part of muffin-tin
        do irc=1,nrci
          call lopzflm(lmaxinr,wfmt1(:,irc,jst),lmmaxvr,zlflm)
          t1=socfr(irc,ias)
          wfmt4(1:lmi,irc,1)=wfmt4(1:lmi,irc,1)+t1*zlflm(1:lmi,3)
          wfmt4(1:lmi,irc,2)=wfmt4(1:lmi,irc,2)-t1*zlflm(1:lmi,3)
          if (ncmag) wfmt4(1:lmi,irc,3)=wfmt4(1:lmi,irc,3)+t1*(zlflm(1:lmi,1) &
           +cmplx(aimag(zlflm(1:lmi,2)),-dble(zlflm(1:lmi,2)),8))
        end do
! outer part of muffin-tin
        do irc=iro,nrc
          call lopzflm(lmaxvr,wfmt1(:,irc,jst),lmmaxvr,zlflm)
          t1=socfr(irc,ias)
          wfmt4(:,irc,1)=wfmt4(:,irc,1)+t1*zlflm(:,3)
          wfmt4(:,irc,2)=wfmt4(:,irc,2)-t1*zlflm(:,3)
          if (ncmag) wfmt4(:,irc,3)=wfmt4(:,irc,3)+t1*(zlflm(:,1) &
           +cmplx(aimag(zlflm(:,2)),-dble(zlflm(:,2)),8))
        end do
      end if
    else
      wfmt4(:,:,:)=0.d0
    end if
! apply muffin-tin potential matrix if required
    if (tvmatmt) then
      do l=0,lmaxdm
        if (tvmmt(l,ias)) then
          nm=2*l+1
          lm=idxlm(l,-l)
          do k=1,nsc
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
            call zgemm('N','N',nm,nrc,nm,zone,vmatmt(lm,ispn,lm,jspn,ias),ld, &
             wfmt1(lm,1,jst),lmmaxvr,zone,wfmt4(lm,1,k),lmmaxvr)
          end do
        end if
      end do
    end if
! apply vector potential if required
    if (afieldpol) then
      call gradzfmt(nrc,nrci,rcmt(:,is),wfmt1(:,:,jst),nrcmtmax,gwfmt)
      do irc=1,nrc
        wfmt3(:,irc)=afieldc(1)*gwfmt(:,irc,1) &
                    +afieldc(2)*gwfmt(:,irc,2) &
                    +afieldc(3)*gwfmt(:,irc,3)
        wfmt3(:,irc)=ca*cmplx(aimag(wfmt3(:,irc)),-dble(wfmt3(:,irc)),8)
      end do
      do k=1,nsd
        wfmt4(:,1:nrc,k)=wfmt4(:,1:nrc,k)+wfmt3(:,1:nrc)
      end do
    end if
! second-variational Hamiltonian matrix
    do ist=1,nstfv
      do k=1,nsc
        if (k.eq.1) then
          i=ist
          j=jst
        else if (k.eq.2) then
          i=ist+nstfv
          j=jst+nstfv
        else
          i=ist
          j=jst+nstfv
        end if
        if (i.le.j) then
!$OMP ATOMIC
          evecsv(i,j)=evecsv(i,j)+zfmtinp(.true.,nrc,nrci,rcmt(:,is), &
           wfmt1(:,:,ist),wfmt4(:,:,k))
        end if
      end do
    end do
! end loop over states
  end do
  deallocate(wfmt1,wfmt2,wfmt3,wfmt4)
  if (afieldpol) deallocate(gwfmt)
! end loop over atoms
end do
!$OMP END DO
!$OMP END PARALLEL
!---------------------------!
!     interstitial part     !
!---------------------------!
if (spinpol) then
! begin loop over states
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,wfir2,z) &
!$OMP PRIVATE(igp,ifg,t1,z1,i,j,k,ist)
!$OMP DO
  do jst=1,nstfv
    allocate(wfir1(ngtot),wfir2(ngtot),z(ngkmax,nsc))
    wfir1(:)=0.d0
    do igp=1,ngp
      ifg=igfft(igpig(igp))
      wfir1(ifg)=evecfv(igp,jst)
    end do
! Fourier transform wavefunction to real-space
    call zfftifc(3,ngridg,1,wfir1)
! multiply with magnetic field and transform to G-space
    wfir2(:)=bsir(:,ndmag)*wfir1(:)
    call zfftifc(3,ngridg,-1,wfir2)
    do igp=1,ngp
      ifg=igfft(igpig(igp))
      z(igp,1)=wfir2(ifg)
      z(igp,2)=-wfir2(ifg)
    end do
    if (ncmag) then
      wfir2(:)=cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:)
      call zfftifc(3,ngridg,-1,wfir2)
      do igp=1,ngp
        ifg=igfft(igpig(igp))
        z(igp,3)=wfir2(ifg)
      end do
    end if
! apply vector potential if required
    if (afieldpol) then
! multiply wavefunction with characteristic function and transform to G-space
      wfir1(:)=wfir1(:)*cfunir(:)
      call zfftifc(3,ngridg,-1,wfir1)
      do igp=1,ngp
        ifg=igfft(igpig(igp))
        t1=ca*dot_product(afieldc(:),vgpc(:,igp))
        z1=t1*wfir1(ifg)
        do k=1,nsd
          z(igp,k)=z(igp,k)+z1
        end do
      end do
    end if
! add to Hamiltonian matrix
    do ist=1,nstfv
      do k=1,nsc
        if (k.eq.1) then
          i=ist
          j=jst
        else if (k.eq.2) then
          i=ist+nstfv
          j=jst+nstfv
        else
          i=ist
          j=jst+nstfv
        end if
        if (i.le.j) then
!$OMP ATOMIC
          evecsv(i,j)=evecsv(i,j)+zdotc(ngp,evecfv(:,ist),1,z(:,k),1)
        end if
      end do
    end do
    deallocate(wfir1,wfir2,z)
! end loop over states
  end do
!$OMP END DO
!$OMP END PARALLEL
end if
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist)
  end do
end do
! diagonalise second-variational Hamiltonian
allocate(rwork(3*nstsv))
lwork=2*nstsv
allocate(work(lwork))
if (ndmag.eq.1) then
! collinear: block diagonalise H
  call zheev('V','U',nstfv,evecsv,nstsv,evalsvp,work,lwork,rwork,info)
  if (info.ne.0) goto 20
  i=nstfv+1
  call zheev('V','U',nstfv,evecsv(i,i),nstsv,evalsvp(i),work,lwork,rwork,info)
  if (info.ne.0) goto 20
  do i=1,nstfv
    do j=1,nstfv
      evecsv(i,j+nstfv)=0.d0
      evecsv(i+nstfv,j)=0.d0
    end do
  end do
else
! non-collinear or spin-unpolarised: full diagonalisation
  call zheev('V','U',nstsv,evecsv,nstsv,evalsvp,work,lwork,rwork,info)
  if (info.ne.0) goto 20
end if
deallocate(rwork,work)
call timesec(ts1)
!$OMP ATOMIC
timesv=timesv+ts1-ts0
return
20 continue
write(*,*)
write(*,'("Error(eveqnsv): diagonalisation of the second-variational &
 &Hamiltonian failed")')
write(*,'(" ZHEEV returned INFO = ",I8)') info
write(*,*)
stop
end subroutine
