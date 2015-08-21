
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlxbsek(ik2,jlgr)
use modmain
implicit none
! arguments
integer, intent(in) :: ik2
real(8), intent(in) :: jlgr(0:lmaxvr+npsden+1,ngvec,nspecies)
! local variables
integer i1,i2,j1,j2
integer a1,a2,b1,b2
integer ik1,is,ia,ias
integer ist1,ist2,jst1,jst2
complex(8) zrho0,zt1
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:)
complex(8), allocatable :: wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
complex(8), allocatable :: zfmt(:,:)
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot))
allocate(zvclir(ngrtot))
allocate(zfmt(lmmaxvr,nrcmtmax))
! calculate the wavefunctions for all states for k-point ik2
call genwfsvp(.false.,.false.,vkl(:,ik2),wfmt2,ngrtot,wfir2)
do i2=1,nvbse
  ist2=istbse(i2,ik2)
  do j2=1,ncbse
    jst2=jstbse(j2,ik2)
    a2=ijkbse(i2,j2,ik2)
! calculate the complex overlap density
    call genzrho(.true.,wfmt2(:,:,:,:,ist2),wfmt2(:,:,:,:,jst2), &
     wfir2(:,:,ist2),wfir2(:,:,jst2),zrhomt,zrhoir)
! compute the potential and G = 0 coefficient of the density
    call genzvclmt(nrcmt,nrcmtmax,rcmt,nrcmtmax,zrhomt,zvclmt)
    call zpotcoul(nrcmt,nrcmtmax,rcmt,1,gc,jlgr,ylmg,sfacg,zrhoir,nrcmtmax, &
     zvclmt,zvclir,zrho0)
! convert potential to spherical coordinates
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        zfmt(:,1:nrcmt(is))=zvclmt(:,1:nrcmt(is),ias)
        call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtvr,lmmaxvr, &
         zfmt,lmmaxvr,zzero,zvclmt(:,:,ias),lmmaxvr)
      end do
    end do
! start loop over ik1
    do ik1=1,nkptnr
      if (ik1.eq.ik2) then
        wfmt1(:,:,:,:,:)=wfmt2(:,:,:,:,:)
        wfir1(:,:,:)=wfir2(:,:,:)
      else
        call genwfsvp(.false.,.false.,vkl(:,ik1),wfmt1,ngrtot,wfir1)
      end if
      do i1=1,nvbse
        ist1=istbse(i1,ik1)
        do j1=1,ncbse
          jst1=jstbse(j1,ik1)
          a1=ijkbse(i1,j1,ik1)
! calculate the complex overlap density
          call genzrho(.false.,wfmt1(:,:,:,:,ist1),wfmt1(:,:,:,:,jst1), &
           wfir1(:,:,ist1),wfir1(:,:,jst1),zrhomt,zrhoir)
! compute the matrix element
          zt1=2.d0*wkptnr*zfinp(.false.,zrhomt,zvclmt,zrhoir,zvclir)
          hmlbse(a1,a2)=hmlbse(a1,a2)+zt1
! compute off-diagonal blocks if required
          if (bsefull) then
            b1=a1+nbbse
            b2=a2+nbbse
            hmlbse(b1,b2)=hmlbse(b1,b2)-conjg(zt1)
! conjugate the potential
            do is=1,nspecies
              do ia=1,natoms(is)
                ias=idxas(ia,is)
                zvclmt(:,1:nrcmt(is),ias)=conjg(zvclmt(:,1:nrcmt(is),ias))
              end do
            end do
            zvclir(:)=conjg(zvclir(:))
            zt1=2.d0*wkptnr*zfinp(.false.,zrhomt,zvclmt,zrhoir,zvclir)
            hmlbse(a1,b2)=hmlbse(a1,b2)+zt1
            hmlbse(b1,a2)=hmlbse(b1,a2)-conjg(zt1)
          end if
        end do
      end do
! end loop over ik1
    end do
! end loop over j2
  end do
! end loop over i2
end do
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir,zvclmt,zvclir,zfmt)
return
end subroutine
