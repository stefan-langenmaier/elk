
! Copyright (C) 2012 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genspchi0
! !INTERFACE:
subroutine genspchi0(ikp,vqpl,expqmt,chi0)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ikp    : k-point from non-reduced set (in,integer)
!   vqpl   : input q-point in lattice coordinates (in,real(3))
!   expqmt : exp(q.r) function in the muffin-tin in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   chi0   : spin-dependent Kohn-Sham response function in G-space
!            (out,complex(ngrpa,4,ngrpa,4,nwrpa))
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
real(8), intent(in) :: vqpl(3)
complex(8), intent(in) :: expqmt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(inout) :: chi0(ngrpa,4,ngrpa,4,nwrpa)
! local variables
integer ispn,ias,is
integer isym,jkp,jkpq
integer ist,jst,irc,iw
integer ig,jg,a,b,i,j
real(8) vkql(3),eij,t1
complex(8) zv(4),zt1,zt2
! allocatable arrays
complex(8), allocatable :: wfmtp(:,:,:,:,:),wfmtpq(:,:,:,:,:)
complex(8), allocatable :: wfirp(:,:,:),wfirpq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:)
complex(8), allocatable :: zw(:),zg(:,:)
! external functions
complex(8) zfinp
external zfinp
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
allocate(wfmtp(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirp(ngrtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,vkl(:,ikp),wfmtp,ngrtot,wfirp)
allocate(wfmtpq(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfirpq(ngrtot,nspinor,nstsv))
call genwfsvp(.false.,.false.,vkql,wfmtpq,ngrtot,wfirpq)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zw,zg) &
!$OMP PRIVATE(ispn,ias,is,irc,jst) &
!$OMP PRIVATE(t1,eij,iw,i,j,a,b) &
!$OMP PRIVATE(ig,jg,zv,zt1,zt2)
!$OMP DO
do ist=1,nstsv
  allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot),zrhoir(ngrtot))
  allocate(zw(nwrpa),zg(ngrpa,4))
! multiply wfmtp with exp(iq.r)
  do ispn=1,nspinor
    do ias=1,natmtot
      is=idxis(ias)
      do irc=1,nrcmt(is)
        wfmtp(:,irc,ias,ispn,ist)=wfmtp(:,irc,ias,ispn,ist)*expqmt(:,irc,ias)
      end do
    end do
  end do
  do jst=1,nstsv
    t1=(wkptnr/omega)*(occsv(ist,jkp)-occsv(jst,jkpq))
    if (abs(t1).gt.1.d-8) then
      eij=evalsv(ist,jkp)-evalsv(jst,jkpq)
! scissor operator
      if (abs(scissor).gt.0.d0) then
        if (eij.gt.0.d0) then
          eij=eij+scissor
        else
          eij=eij-scissor
        end if
      end if
! frequency-dependent part in response function formula for all frequencies
      do iw=1,nwrpa
        zw(iw)=t1/(eij+wrpa(iw))
      end do
! compute the complex density in G+q-space
      i=0
      do a=1,2
        do b=1,2
          i=i+1
          call genzrho(.false.,.false.,wfmtp(:,:,:,b,ist),wfmtpq(:,:,:,a,jst), &
           wfirp(:,b,ist),wfirpq(:,a,jst),zrhomt,zrhoir)
          do ig=1,ngrpa
            zg(ig,i)=zfinp(.false.,zrhomt,expgmt(:,:,:,ig),zrhoir,expgir(:,ig))
          end do
        end do
      end do
!$OMP CRITICAL
      do ig=1,ngrpa
! conjugate transpose of matrix
        zv(1)=conjg(zg(ig,1))
        zv(2)=conjg(zg(ig,2))
        zv(3)=conjg(zg(ig,3))
        zv(4)=conjg(zg(ig,4))
        do j=1,4
          do jg=1,ngrpa
            zt1=zg(jg,j)
            do i=1,4
              zt2=zv(i)*zt1
              chi0(ig,i,jg,j,:)=chi0(ig,i,jg,j,:)+zt2*zw(:)
            end do
          end do
        end do
      end do
!$OMP END CRITICAL
    end if
! end loop over jst
  end do
  deallocate(zrhomt,zrhoir,zw,zg)
! end loop over ist
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(wfmtp,wfmtpq,wfirp,wfirpq)
return
end subroutine
!EOC

