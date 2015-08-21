
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine alpha2f
use modmain
use modphonon
use modtest
implicit none
! local variables
integer nb,ik,iq,i,j
integer i1,i2,i3,iw
integer lwork,info
real(8) wmin,wmax,wd,dw
real(8) wlog,wrms,lambda,tc
real(8) v(3),t1
! allocatable arrays
real(8), allocatable :: wq(:,:),wp(:),gq(:,:),a2fp(:)
real(8), allocatable :: w(:),a2f(:),rwork(:)
complex(8), allocatable :: dynq(:,:,:),dynr(:,:,:)
complex(8), allocatable :: dynp(:,:),ev(:,:),b(:,:)
complex(8), allocatable :: a2fmq(:,:,:),a2fmr(:,:,:),a2fmp(:,:)
complex(8), allocatable :: work(:)
! initialise universal variables
call init0
call init1
call init2
nb=3*natmtot
allocate(wq(nb,nqpt),wp(nb),gq(nb,nqpt),a2fp(nb))
allocate(w(nwplot),a2f(nwplot),rwork(3*nb))
allocate(dynq(nb,nb,nqpt))
allocate(dynr(nb,nb,nqptnr))
allocate(dynp(nb,nb))
allocate(ev(nb,nb),b(nb,nb))
allocate(a2fmq(nb,nb,nqpt))
allocate(a2fmr(nb,nb,nqptnr))
allocate(a2fmp(nb,nb))
lwork=2*nb
allocate(work(lwork))
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! compute the density of states at the Fermi energy
call occupy
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! Fourier transform the dynamical matrices to real-space
call dynqtor(dynq,dynr)
! read in the phonon linewidths for each q-point
call readgamma(gq)
! loop over phonon q-points
do iq=1,nqpt
! diagonalise the dynamical matrix
  call dyndiag(dynq(:,:,iq),wq(:,iq),ev)
! construct a complex matrix from the phonon eigenvectors such that its
! eigenvalues squared are the phonon linewidths divided by the frequency
  do i=1,nb
    if (wq(i,iq).gt.1.d-8) then
      t1=sqrt(abs(gq(i,iq)/wq(i,iq)))
    else
      t1=0.d0
    end if
    do j=1,nb
      b(i,j)=t1*conjg(ev(j,i))
    end do
  end do
  call zgemm('N','N',nb,nb,nb,zone,ev,nb,b,nb,zzero,a2fmq(:,:,iq),nb)
end do
! Fourier transform the matrices to real-space
call dynqtor(a2fmq,a2fmr)
! find the minimum and maximum frequencies
wmin=0.d0
wmax=0.d0
do iq=1,nqpt
  wmin=min(wmin,wq(1,iq))
  wmax=max(wmax,wq(nb,iq))
end do
wmax=wmax+(wmax-wmin)*0.1d0
wmin=wmin-(wmax-wmin)*0.1d0
wd=wmax-wmin
if (wd.lt.1.d-8) wd=1.d0
dw=wd/dble(nwplot)
! generate energy grid
do iw=1,nwplot
  w(iw)=dw*dble(iw-1)+wmin
end do
a2f(:)=0.d0
do i1=0,ngrkf-1
  v(1)=dble(i1)/dble(ngrkf)
  do i2=0,ngrkf-1
    v(2)=dble(i2)/dble(ngrkf)
    do i3=0,ngrkf-1
      v(3)=dble(i3)/dble(ngrkf)
! compute the dynamical matrix at this particular q-point
      call dynrtoq(v,dynr,dynp)
! find the phonon frequencies
      call dyndiag(dynp,wp,ev)
! compute the alpha^2F matrix at this particular q-point
      call dynrtoq(v,a2fmr,a2fmp)
! diagonlise the alpha^2F matrix
      call zheev('N','U',nb,a2fmp,nb,a2fp,work,lwork,rwork,info)
! square the eigenvalues to recover the linewidths divided by the frequency
      a2fp(:)=a2fp(:)**2
      do i=1,nb
        t1=(wp(i)-wmin)/dw+1.d0
        iw=nint(t1)
        if ((iw.ge.1).and.(iw.le.nwplot)) then
          a2f(iw)=a2f(iw)+a2fp(i)
        end if
      end do
    end do
  end do
end do
t1=twopi*(fermidos/2.d0)*dw*dble(ngrkf)**3
if (t1.gt.1.d-8) then
  t1=1.d0/t1
else
  t1=0.d0
end if
a2f(:)=t1*a2f(:)
! smooth Eliashberg function if required
if (nswplot.gt.0) call fsmooth(nswplot,nwplot,1,a2f)
! write Eliashberg function to file
open(50,file='ALPHA2F.OUT',action='WRITE',form='FORMATTED')
do iw=1,nwplot
  write(50,'(2G18.10)') w(iw),a2f(iw)
end do
close(50)
write(*,*)
write(*,'("Info(alpha2f):")')
write(*,'(" Eliashberg function written to ALPHA2F.OUT")')
! compute lambda, logarithmic average frequency, RMS average frequency and
! McMillan-Allen-Dynes superconducting critical temperature
call mcmillan(w,a2f,lambda,wlog,wrms,tc)
open(50,file='MCMILLAN.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("Electron-phonon coupling constant, lambda : ",G18.10)') lambda
write(50,*)
write(50,'("Logarithmic average frequency : ",G18.10)') wlog
write(50,*)
write(50,'("RMS average frequency : ",G18.10)') wrms
write(50,*)
write(50,'("Coulomb pseudopotential, mu* : ",G18.10)') mustar
write(50,*)
write(50,'("McMillan-Allen-Dynes superconducting critical temperature")')
write(50,'(" [Eq. 34, Phys. Rev. B 12, 905 (1975)] (kelvin) : ",G18.10)') tc
write(50,*)
close(50)
write(*,*)
write(*,'("Info(alpha2f):")')
write(*,'(" Electron-phonon coupling constant, lambda;")')
write(*,'(" logarithmic and RMS average frequencies;")')
write(*,'(" and McMillan-Allen-Dynes superconducting critical temperature")')
write(*,'(" written to MCMILLAN.OUT")')
! write lambda to test file
call writetest(251,'Electron-phonon coupling constant, lambda',tol=5.d-2, &
 rv=lambda)
deallocate(wq,wp,gq,a2fp,w,a2f)
deallocate(rwork,dynq,dynr,dynp,ev,b)
deallocate(a2fmq,a2fmr,a2fmp,work)
return
end subroutine

