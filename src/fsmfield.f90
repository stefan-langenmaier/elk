
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
integer is,ia,ias
integer idm,ir,irc
real(8) v(3),t1,t2,t3
if ((.not.spinpol).or.(fixspin.eq.0)) return
v(:)=0.d0
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
      t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
      t2=sqrt(momtot(1)**2+momtot(2)**2+momtot(3)**2)
    else
      t1=abs(v(1))
      t2=abs(momtot(1))
    end if
    if (t1.gt.1.d-8) t1=t2/t1
    v(1:ndmag)=t1*v(1:ndmag)
  end if
  do idm=1,ndmag
    bfsmc(idm)=bfsmc(idm)+taufsm*(momtot(idm)-v(idm))
  end do
  do idm=1,ndmag
    t3=bfsmc(idm)
    do ias=1,natmtot
      is=idxis(ias)
      do irc=1,nrcmt(is)
        bsmt(:,irc,ias,idm)=bsmt(:,irc,ias,idm)+t3
      end do
    end do
    do ir=1,ngtot
      bsir(ir,idm)=bsir(ir,idm)+t3*cfunir(ir)
    end do
  end do
end if
if ((abs(fixspin).eq.2).or.(abs(fixspin).eq.3)) then
! determine the muffin-tin fields for fixed local moments
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! if any component is >= 1000 then do not fix the moment
      t1=sum(abs(mommtfix(:,ia,is)))
      if (t1.ge.1000.d0) goto 10
      if (ncmag) then
        v(:)=mommtfix(:,ia,is)
      else
        v(1)=mommtfix(3,ia,is)
      end if
! case where only the direction is fixed
      if (fixspin.lt.0) then
        if (ncmag) then
          t1=sqrt(v(1)**2+v(2)**2+v(3)**2)
          t2=sqrt(mommt(1,ias)**2+mommt(2,ias)**2+mommt(3,ias)**2)
        else
          t1=abs(v(1))
          t2=abs(mommt(1,ias))
        end if
        if (t1.gt.1.d-8) t1=t2/t1
        v(1:ndmag)=t1*v(1:ndmag)
      end if
      do idm=1,ndmag
        bfsmcmt(idm,ias)=bfsmcmt(idm,ias)+taufsm*(mommt(idm,ias)-v(idm))
        t3=bfsmcmt(idm,ias)
        do irc=1,nrcmt(is)
          bsmt(:,irc,ias,idm)=bsmt(:,irc,ias,idm)+t3
        end do
      end do
10 continue
    end do
  end do
end if
return
end subroutine
!EOC

