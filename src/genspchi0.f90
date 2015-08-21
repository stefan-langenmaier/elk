
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genspchi0
! !INTERFACE:
subroutine genspchi0(ikp,scsr,vqpl,gqc,ylmgq,sfacgq,chi0)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp    : k-point from non-reduced set (in,integer)
!   scsr   : scissor correction (in,real)
!   vqpl   : input q-point in lattice coordinates (in,real(3))
!   gqc    : length of G+q-vectors (in,real(ngrf))
!   ylmgq  : spherical harmonics of the G+q-vectors (in,complex(lmmaxvr,ngrf))
!   sfacgq : structure factors of G+q-vectors (in,complex(ngrf,natmtot))
!   chi0   : spin-dependent Kohn-Sham response function in G-space
!            (out,complex(ngrf,4,ngrf,4,nwrf))
! !DESCRIPTION:
!   Computes the spin-dependent Kohn-Sham response function:
!   \begin{align*}
!    \chi_{\alpha\beta,\alpha'\beta'}({\bf r},{\bf r}',\omega)
!    & \equiv\frac{\delta\rho_{\alpha\beta}({\bf r},\omega)}
!    {\delta v_{\alpha'\beta'}({\bf r}',\omega)} \\
!    & =\frac{1}{N_k}\sum_{i{\bf k},j{\bf k}'}(f_{i{\bf k}}-f_{j{\bf k}'})
!    \frac{\langle i{\bf k}|\hat{\rho}_{\beta\alpha}({\bf r})|j{\bf k}'\rangle
!    \langle j{\bf k}'|\hat{\rho}_{\alpha'\beta'}({\bf r}')|i{\bf k}\rangle}
!    {w+(\varepsilon_{i{\bf k}}-\varepsilon_{j{\bf k}'})+i\eta},
!   \end{align*}
!   where $\alpha$ and $\beta$ are spin-coordinates, $N_k$ is the number of
!   $k$-points, $f_{i{\bf k}}$ are the occupancies, $v$ is the Kohn-Sham
!   potential and $\hat{\rho}$ is the spin-density operator. With translational
!   symmetry in mind, we adopt the following convention for its Fourier
!   transform:
!   $$ \chi_{\alpha\beta,\alpha'\beta'}({\bf G},{\bf G}',{\bf q},\omega)=
!    \frac{1}{\Omega}\int d^3r\,d^3r'\,e^{-i({\bf G}+{\bf q})\cdot{\bf r}}
!    e^{i({\bf G}'+{\bf q})\cdot{\bf r}'}
!    \chi_{\alpha\beta,\alpha'\beta'}({\bf r},{\bf r}',\omega). $$
!   Let
!   $$ Z_{i{\bf k},j{\bf k}+{\bf q}}^{\alpha\beta}({\bf G})\equiv
!    \int d^3r\,e^{i({\bf G}+{\bf q})\cdot{\bf r}}
!    \varphi_{j{\bf k}+{\bf q},\alpha}^*({\bf r})
!    \varphi_{i{\bf k},\beta}({\bf r}) $$
!   then the response function in $G$-space can be written
!   $$ \chi_{\alpha\beta,\alpha'\beta'}({\bf G},{\bf G}',{\bf q},\omega)=
!    \frac{1}{N_k\Omega}\sum_{i{\bf k},j{\bf k}+{\bf q}}
!    (f_{i{\bf k}}-f_{j{\bf k}})
!    \frac{\left[Z_{i{\bf k},j{\bf k}+{\bf q}}^{\alpha\beta}({\bf G})\right]^*
!    Z_{i{\bf k},j{\bf k}+{\bf q}}^{\alpha'\beta'}({\bf G}')}
!    {w+(\varepsilon_{i{\bf k}}-\varepsilon_{j{\bf k}+{\bf q}})+i\eta}. $$
!
! !REVISION HISTORY:
!   Created March 2012 (SS and JKD)
!EOP
!BOC
implicit none
! local variables
integer, intent(in) :: ikp
real(8), intent(in) :: scsr
real(8), intent(in) :: vqpl(3)
real(8), intent(in) :: gqc(ngrf)
complex(8), intent(in) :: ylmgq(lmmaxvr,ngrf)
complex(8), intent(in) :: sfacgq(ngrf,natmtot)
complex(8), intent(inout) :: chi0(nwrf,ngrf,4,ngrf,4)
! local variables
logical tz(4)
integer isym,jkp,jkpq,iw
integer ig,jg,ist,jst,a,b,i,j
real(8) vkql(3),eij,t1
complex(8) z1,z2
! allocatable arrays
complex(8), allocatable :: wfmt(:,:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:,:),wfirq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zw(:),zg(:,:)
if (.not.spinpol) then
  write(*,*)
  write(*,'("Error(genspchi0): spin-unpolarised calculation")')
  write(*,*)
  stop
end if
! k+q-vector in lattice coordinates
vkql(:)=vkl(:,ikp)+vqpl(:)
! equivalent reduced k-points for k and k+q
call findkpt(vkl(:,ikp),isym,jkp)
call findkpt(vkql,isym,jkpq)
! generate the wavefunctions for all states at k and k+q
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,.false.,vkl(:,ikp),wfmt,ngtot,wfir)
allocate(wfmtq(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirq(ngtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,.false.,vkql,wfmtq,ngtot,wfirq)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zw,zg) &
!$OMP PRIVATE(jst,t1,eij,iw,i,j) &
!$OMP PRIVATE(a,b,tz,ig,jg,z1,z2)
!$OMP DO
do ist=1,nstsv
  if (abs(evalsv(ist,jkp)-efermi).gt.emaxrf) cycle
  allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngtot))
  allocate(zw(nwrf),zg(ngrf,4))
  do jst=1,nstsv
    if (abs(evalsv(jst,jkpq)-efermi).gt.emaxrf) cycle
    t1=wkptnr*omega*(occsv(ist,jkp)-occsv(jst,jkpq))
    if (abs(t1).lt.1.d-8) cycle
    eij=evalsv(ist,jkp)-evalsv(jst,jkpq)
! scissor operator
    if (abs(scsr).gt.1.d-8) then
      if (eij.gt.0.d0) then
        eij=eij+scsr
      else
        eij=eij-scsr
      end if
    end if
! frequency-dependent part in response function formula for all frequencies
    do iw=1,nwrf
      zw(iw)=t1/(eij+wrf(iw))
    end do
! compute the complex density in G+q-space
    i=0
    do a=1,2
      do b=1,2
        i=i+1
! find which contributions are zero for collinear case
        tz(i)=.false.
        if (.not.ncmag) then
          if (((a.eq.1).and.(ist.gt.nstfv)).or. &
              ((a.eq.2).and.(ist.le.nstfv)).or. &
              ((b.eq.1).and.(jst.gt.nstfv)).or. &
              ((b.eq.2).and.(jst.le.nstfv))) then
            tz(i)=.true.
            cycle
          end if
        end if
        call genzrho(.true.,.false.,wfmt(:,:,:,a,ist),wfmtq(:,:,:,b,jst), &
         wfir(:,a,ist),wfirq(:,b,jst),zrhomt,zrhoir)
        call zftzf(ngrf,gqc,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zg(:,i))
      end do
    end do
!$OMP CRITICAL
    do j=1,4
      if (tz(j)) cycle
      do jg=1,ngrf
        z1=conjg(zg(jg,j))
        do i=1,4
          if (tz(i)) cycle
          do ig=1,ngrf
            z2=zg(ig,i)*z1
            call zaxpy(nwrf,z2,zw,1,chi0(:,ig,i,jg,j),1)
          end do
        end do
      end do
    end do
!$OMP END CRITICAL
! end loop over jst
  end do
  deallocate(zrhomt,zrhoir,zw,zg)
! end loop over ist
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(wfmt,wfmtq,wfir,wfirq)
return
end subroutine
!EOC

