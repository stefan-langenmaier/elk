
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: forcek
! !INTERFACE:
subroutine forcek(ik)
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the {\bf k}-dependent contribution to the incomplete basis set
!   (IBS) force. See the calling routine {\tt force} for a full description.
!
! !REVISION HISTORY:
!   Created June 2006 (JKD)
!   Updated for spin-spiral case, May 2007 (Francesco Cricchio and JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ispn,jspn,jspn0,jspn1
integer is,ia,ias,ist,jst
integer n2,i,j,k,l,iv(3),ig
real(8) sum,t1
complex(8) zt1,zt2
! allocatable arrays
integer, allocatable :: ijg(:)
real(8), allocatable :: dp(:)
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: h(:),o(:),dlh(:),dlo(:)
complex(8), allocatable :: vh(:),vo(:)
complex(8), allocatable :: ffv(:,:),y(:)
! external functions
complex(8) zdotc
external zdotc
n2=nmat(1,ik)**2
if (spinsprl) n2=max(n2,nmat(2,ik)**2)
! allocate local arrays
allocate(ijg(n2))
allocate(dp(n2))
allocate(evalfv(nstfv,nspnfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(h(n2),o(n2),dlh(n2),dlo(n2))
allocate(vh(nmatmax),vo(nmatmax))
allocate(ffv(nstfv,nstfv),y(nstfv))
! get the eigenvalues/vectors and occupancies from file
call getevalfv(vkl(:,ik),evalfv)
call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(vkl(:,ik),evecsv)
call getoccsv(vkl(:,ik),occsv(:,ik))
! loop over first-variational spin components
do ispn=1,nspnfv
  if (spinsprl) then
    jspn0=ispn; jspn1=ispn
  else
    jspn0=1; jspn1=nspinor
  end if
! find the matching coefficients
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm)
  do j=1,ngk(ispn,ik)
    k=(j-1)*nmat(ispn,ik)
    do i=1,j
      k=k+1
      iv(:)=ivg(:,igkig(i,ispn,ik))-ivg(:,igkig(j,ispn,ik))
      iv(:)=modulo(iv(:)-intgv(:,1),ngrid(:))+intgv(:,1)
      ijg(k)=ivgig(iv(1),iv(2),iv(3))
      dp(k)=0.5d0*dot_product(vgkc(:,i,ispn,ik),vgkc(:,j,ispn,ik))
    end do
  end do
! loop over species and atoms
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      h(1:n2)=0.d0
      call hmlaa(ias,ngk(ispn,ik),apwalm,h)
      call hmlalo(ias,ngk(ispn,ik),apwalm,h)
      o(1:n2)=0.d0
      call olpaa(ias,ngk(ispn,ik),apwalm,o)
      call olpalo(ias,ngk(ispn,ik),apwalm,o)
! loop over Cartesian directions
      do l=1,3
! APW-APW contribution
        do j=1,ngk(ispn,ik)
          k=(j-1)*nmat(ispn,ik)
          do i=1,j
            k=k+1
            ig=ijg(k)
            t1=vgc(l,ig)
            zt1=-ffacg(ig,is)*conjg(sfacg(ig,ias))
            dlh(k)=(dp(k)*zt1+h(k))*t1
            dlo(k)=(zt1+o(k))*t1
          end do
        end do
! APW-local-orbital contribution
        do j=ngk(ispn,ik)+1,nmat(ispn,ik)
          k=(j-1)*nmat(ispn,ik)
          do i=1,ngk(ispn,ik)
            k=k+1
            t1=vgkc(l,i,ispn,ik)
            dlh(k)=h(k)*t1
            dlo(k)=o(k)*t1
          end do
          do i=ngk(ispn,ik)+1,j
            k=k+1
            dlh(k)=0.d0
            dlo(k)=0.d0
          end do
        end do
! multiply by i
        do k=1,n2
          dlh(k)=cmplx(-aimag(dlh(k)),dble(dlh(k)),8)
          dlo(k)=cmplx(-aimag(dlo(k)),dble(dlo(k)),8)
        end do
! compute the force matrix elements in the first-variational basis
        do jst=1,nstfv
          call zhemv('U',nmat(ispn,ik),zone,dlh,nmat(ispn,ik), &
           evecfv(:,jst,ispn),1,zzero,vh,1)
          call zhemv('U',nmat(ispn,ik),zone,dlo,nmat(ispn,ik), &
           evecfv(:,jst,ispn),1,zzero,vo,1)
          t1=evalfv(jst,ispn)
          do ist=1,nstfv
            zt1=zdotc(nmat(ispn,ik),evecfv(:,ist,ispn),1,vh,1)
            zt2=zdotc(nmat(ispn,ik),evecfv(:,ist,ispn),1,vo,1)
            ffv(ist,jst)=zt1-t1*zt2
          end do
        end do
! compute the force using the second-variational coefficients if required
        sum=0.d0
        if (tevecsv) then
! spin-polarised case
          do j=1,nstsv
            t1=occsv(j,ik)
            do jspn=jspn0,jspn1
              i=(jspn-1)*nstfv+1
              call zgemv('N',nstfv,nstfv,zone,ffv,nstfv,evecsv(i,j),1,zzero,y,1)
              zt1=zdotc(nstfv,evecsv(i,j),1,y,1)
              sum=sum+t1*dble(zt1)
            end do
          end do
        else
! spin-unpolarised case
          do j=1,nstsv
            sum=sum+occsv(j,ik)*dble(ffv(j,j))
          end do
        end if
!$OMP CRITICAL
        forceibs(l,ias)=forceibs(l,ias)+wkpt(ik)*sum
!$OMP END CRITICAL
! end loop over Cartesian components
      end do
! end loop over atoms and species
    end do
  end do
! end loop over first-variational spins
end do
deallocate(ijg,dp,evalfv,apwalm,evecfv,evecsv)
deallocate(h,o,dlh,dlo,vh,vo,ffv,y)
return
end subroutine
!EOC

