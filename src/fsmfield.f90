
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: fsmfield
! !INTERFACE:
subroutine fsmfield
! !USES:
use modmain
! !DESCRIPTION:
!   Updates the effective magnetic field, ${\bf B}_{\rm FSM}$, required for
!   fixing the spin moment to a given value, $\boldsymbol{\mu}_{\rm FSM}$. This
!   is done by adding a vector to the field which is proportional to the
!   difference between the moment calculated in the $i$th self-consistent loop
!   and the required moment:
!   $$ {\bf B}_{\rm FSM}^{i+1}={\bf B}_{\rm FSM}^i+\lambda\left(
!    \boldsymbol{\mu}^i-\boldsymbol{\mu}_{\rm FSM}\right), $$
!   where $\lambda$ is a scaling factor.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir,idm
real(8) v(3),t1,t2,t3
if ((.not.spinpol).or.(fixspin.eq.0)) return
v(:)=0.d0
t1=1.d0/y00
! determine the global effective field
if ((abs(fixspin).eq.1).or.(abs(fixspin).eq.3)) then
  if (ncmag) then
    v(:)=momfix(:)
  else
    v(1)=momfix(3)
  end if
! case where only the direction is fixed
  if (fixspin.lt.0) then
    if (ncmag) then
      t2=sqrt(v(1)**2+v(2)**2+v(3)**2)
      t3=sqrt(momtot(1)**2+momtot(2)**2+momtot(3)**2)
    else
      t2=abs(v(1))
      t3=abs(momtot(1))
    end if
    if (t2.gt.1.d-8) t2=t3/t2
    v(1:ndmag)=t2*v(1:ndmag)
  end if
  do idm=1,ndmag
    bfsmc(idm)=bfsmc(idm)+taufsm*(momtot(idm)-v(idm))
  end do
  do idm=1,ndmag
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ir=1,nrmt(is)
          bxcmt(1,ir,ias,idm)=bxcmt(1,ir,ias,idm)+t1*bfsmc(idm)
        end do
      end do
    end do
    do ir=1,ngrtot
      bxcir(ir,idm)=bxcir(ir,idm)+bfsmc(idm)
    end do
  end do
end if
if ((abs(fixspin).eq.2).or.(abs(fixspin).eq.3)) then
! determine the muffin-tin fields for fixed local moments
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! if any component is >= 1000 then do not fix the moment
      t2=sum(abs(mommtfix(:,ia,is)))
      if (t2.ge.1000.d0) goto 10
      if (ncmag) then
        v(:)=mommtfix(:,ia,is)
      else
        v(1)=mommtfix(3,ia,is)
      end if
! case where only the direction is fixed
      if (fixspin.lt.0) then
        if (ncmag) then
          t2=sqrt(v(1)**2+v(2)**2+v(3)**2)
          t3=sqrt(mommt(1,ias)**2+mommt(2,ias)**2+mommt(3,ias)**2)
        else
          t2=abs(v(1))
          t3=abs(mommt(1,ias))
        end if
        if (t2.gt.1.d-8) t2=t3/t2
        v(1:ndmag)=t2*v(1:ndmag)
      end if
      do idm=1,ndmag
        bfsmcmt(idm,ia,is)=bfsmcmt(idm,ia,is)+taufsm*(mommt(idm,ias)-v(idm))
        do ir=1,nrmt(is)
          bxcmt(1,ir,ias,idm)=bxcmt(1,ir,ias,idm)+t1*bfsmcmt(idm,ia,is)
        end do
      end do
10 continue
    end do
  end do
end if
return
end subroutine
!EOC

