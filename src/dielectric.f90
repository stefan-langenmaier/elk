
! Copyright (C) 2002-2009 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dielectric
! !INTERFACE:
subroutine dielectric
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Computes the dielectric tensor, optical conductivity and plasma frequency.
!   The formulae are taken from {\it Physica Scripta} {\bf T109}, 170 (2004).
!
! !REVISION HISTORY:
!   Created November 2005 (SS and JKD)
!   Added plasma frequency and intraband contribution (S. Lebegue)
!   Complete rewrite, 2008 (JKD)
!   Fixed problem with plasma frequency, 2009 (Marty Blaber and JKD)
!   Parallelised, 2009 (M. Blaber)
!EOP
!BOC
implicit none
! local variables
integer ik,jk,ist,jst
integer iw,i,j,l
real(8) w1,w2,wplas
real(8) eji,x,t1,t2
complex(8) eta,zt1
character(256) fname
! allocatable arrays
real(8), allocatable :: w(:)
real(8), allocatable :: delta(:,:,:)
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: sigma(:)
! external functions
real(8) sdelta
external sdelta
! initialise universal variables
call init0
call init1
! read Fermi energy from file
call readfermi
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! allocate local arrays
allocate(w(nwdos))
if (usegdft) allocate(delta(nstsv,nstsv,nkpt))
allocate(sigma(nwdos))
! compute generalised DFT correction
if (usegdft) then
  call readstate
  call poteff
  call linengy
  call genapwfr
  call genlofr
  do ik=1,nkpt
    call gdft(ik,delta(:,:,ik))
  end do
end if
! generate energy grid (always non-negative)
w1=max(wdos(1),0.d0)
w2=max(wdos(2),w1)
t1=(w2-w1)/dble(nwdos)
do iw=1,nwdos
  w(iw)=w1+t1*dble(iw-1)
end do
! i divided by the complex relaxation time
eta=cmplx(0.d0,swidth,8)
! loop over dielectric tensor components
do l=1,noptcomp
  i=optcomp(1,l)
  j=optcomp(2,l)
  wplas=0.d0
  sigma(:)=0.d0
! parallel loop over non-reduced k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(pmat,jk,ist,jst) &
!$OMP PRIVATE(zt1,eji,t1,x)
!$OMP DO
  do ik=1,nkptnr
    allocate(pmat(3,nstsv,nstsv))
! equivalent reduced k-point
    jk=ikmap(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! read in the momentum matrix elements
    call getpmat(vkl(:,ik),pmat)
! valance states
    do ist=1,nstsv
! conduction states
      do jst=1,nstsv
        zt1=pmat(i,ist,jst)*conjg(pmat(j,ist,jst))
        eji=evalsv(jst,jk)-evalsv(ist,jk)
        if ((evalsv(ist,jk).le.efermi).and.(evalsv(jst,jk).gt.efermi)) then
! generalised DFT correction
          if (usegdft) then
            eji=eji+delta(jst,ist,jk)
          end if
! scissor correction
          if (scissor.ne.0.d0) then
            t1=(eji+scissor)/eji
            zt1=zt1*t1**2
            eji=eji+scissor
          end if
        end if
        if (abs(eji).gt.1.d-8) then
          t1=occsv(ist,jk)*(1.d0-occsv(jst,jk)/occmax)/eji
!$OMP CRITICAL
          sigma(:)=sigma(:)+t1*(zt1/(w(:)-eji+eta)+conjg(zt1)/(w(:)+eji+eta))
!$OMP END CRITICAL
        end if
! add to the plasma frequency
        if (intraband) then
          if (i.eq.j) then
            if (ist.eq.jst) then
              x=(evalsv(ist,jk)-efermi)/swidth
!$OMP CRITICAL
              wplas=wplas+wkptnr*dble(zt1)*sdelta(stype,x)/swidth
!$OMP END CRITICAL
            end if
          end if
        end if
      end do
    end do
    deallocate(pmat)
  end do
!$OMP END DO
!$OMP END PARALLEL
  zt1=zi*wkptnr/omega
  sigma(:)=zt1*sigma(:)
! intraband contribution
  if (intraband) then
    if (i.eq.j) then
      wplas=sqrt(occmax*abs(wplas)*fourpi/omega)
! write the plasma frequency to file
      write(fname,'("PLASMA_",2I1,".OUT")') i,j
      open(60,file=trim(fname),action='WRITE',form='FORMATTED')
      write(60,'(G18.10," : plasma frequency")') wplas
      close(60)
! add the intraband contribution to sigma
      t1=wplas**2/fourpi
      do iw=1,nwdos
        sigma(iw)=sigma(iw)+t1/(swidth-zi*w(iw))
      end do
    end if
  end if
! write the optical conductivity to file
  write(fname,'("SIGMA_",2I1,".OUT")') i,j
  open(60,file=trim(fname),action='WRITE',form='FORMATTED')
  do iw=1,nwdos
    write(60,'(2G18.10)') w(iw),dble(sigma(iw))
  end do
  write(60,'("     ")')
  do iw=1,nwdos
    write(60,'(2G18.10)') w(iw),aimag(sigma(iw))
  end do
  close(60)
! write the dielectric function to file
  write(fname,'("EPSILON_",2I1,".OUT")') i,j
  open(60,file=trim(fname),action='WRITE',form='FORMATTED')
  t1=0.d0
  if (i.eq.j) t1=1.d0
  do iw=1,nwdos
    t2=t1-fourpi*aimag(sigma(iw)/(w(iw)+eta))
    write(60,'(2G18.10)') w(iw),t2
  end do
  write(60,'("     ")')
  do iw=1,nwdos
    t2=fourpi*dble(sigma(iw)/(w(iw)+eta))
    write(60,'(2G18.10)') w(iw),t2
  end do
  close(60)
! write sigma to test file
  call writetest(121,'optical conductivity',nv=nwdos,tol=1.d-2,zva=sigma)
! end loop over tensor components
end do
close(50)
write(*,*)
write(*,'("Info(dielectric):")')
write(*,'(" dielectric tensor written to EPSILON_ij.OUT")')
write(*,'(" optical conductivity written to SIGMA_ij.OUT")')
if (intraband) then
  write(*,'(" plasma frequency written to PLASMA_ij.OUT")')
end if
write(*,'(" for components")')
do l=1,noptcomp
  write(*,'("  i = ",I1,", j = ",I1)') optcomp(1:2,l)
end do
deallocate(w,sigma)
if (usegdft) deallocate(delta)
return
end subroutine
!EOC

