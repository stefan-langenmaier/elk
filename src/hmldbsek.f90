
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmldbsek(ik2)
use modmain
implicit none
! arguments
integer, intent(in) :: ik2
! local variables
integer ik1,ig,jg,n
integer is,ia,ias,irc
integer i1,i2,j1,j2,a1,a2,b1,b2
integer ist1,ist2,jst1,jst2
real(8) vl(3),vc(3)
real(8) vgqc(3),t0,t1
complex(8) zsum
character(256) fname
! allocatable arrays
real(8), allocatable :: gqc(:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: expqmt(:,:,:),expgqmt(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:,:),zrhoir(:,:)
complex(8), allocatable :: zvv(:,:,:),zcc(:,:,:)
complex(8), allocatable :: zvc(:,:,:),zcv(:,:,:)
complex(8), allocatable :: epsinv(:,:,:)
! external functions
complex(8) zfinp
external zfinp
allocate(gqc(ngrpa))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv),wfir2(ngrtot,nspinor,nstsv))
allocate(expqmt(lmmaxvr,nrcmtmax,natmtot))
allocate(expgqmt(lmmaxvr,nrcmtmax,natmtot))
n=max(nvbse,ncbse)
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot,n),zrhoir(ngrtot,n))
allocate(zvv(ngrpa,nvbse,nvbse),zcc(ngrpa,ncbse,ncbse))
allocate(epsinv(ngrpa,ngrpa,nwrpa))
if (bsefull) then
  allocate(zvc(ngrpa,nvbse,ncbse))
  allocate(zcv(ngrpa,ncbse,nvbse))
end if
! generate the wavefunctions at k-point ik2
call genwfsvp(.false.,.false.,vkl(:,ik2),wfmt2,ngrtot,wfir2)
! filename for inverse dielectric function
fname='EPSINV_RPA.OUT'
! begin loop over ik1
do ik1=1,nkptnr
! generate the wavefunctions at k-point ik1
  call genwfsvp(.false.,.false.,vkl(:,ik1),wfmt1,ngrtot,wfir1)
! q-vector in lattice and Cartesian coordinates
  vl(:)=vkl(:,ik1)-vkl(:,ik2)
  vc(:)=vkc(:,ik1)-vkc(:,ik2)
! generate the function exp(iq.r) in the muffin-tins
  call genexpmt(vc,expqmt)
! loop over MBPT G-vectors
  do ig=1,ngrpa
! G+q-vector in Cartesian coordinates
    vgqc(:)=vgc(:,ig)+vc(:)
! length of G+q-vector
    gqc(ig)=sqrt(vgqc(1)**2+vgqc(2)**2+vgqc(3)**2)
! compute the function exp(i(G+q).r) in the muffin-tins
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do irc=1,nrcmt(is)
          expgqmt(:,irc,ias)=expgmt(:,irc,ias,ig)*expqmt(:,irc,ias)
        end do
      end do
    end do
! compute the <v|exp(i(G+q).r)|v'> matrix elements
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,ist2,i2)
!$OMP DO
    do i1=1,nvbse
      ist1=istbse(i1,ik1)
      do i2=1,nvbse
        ist2=istbse(i2,ik2)
! note that the complex conjugate of the density is found because zfinp
! conjugates the first function
        call genzrho(.false.,wfmt2(:,:,:,:,ist2),wfmt1(:,:,:,:,ist1), &
         wfir2(:,:,ist2),wfir1(:,:,ist1),zrhomt(:,:,:,i1),zrhoir(:,i1))
        zvv(ig,i1,i2)=zfinp(.false.,zrhomt(:,:,:,i1),expgqmt,zrhoir(:,i1), &
         expgir(:,ig))
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
! compute the <c|exp(i(G+q).r)|c'> matrix elements
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jst1,jst2,j2)
!$OMP DO
    do j1=1,ncbse
      jst1=jstbse(j1,ik1)
      do j2=1,ncbse
        jst2=jstbse(j2,ik2)
        call genzrho(.false.,wfmt2(:,:,:,:,jst2),wfmt1(:,:,:,:,jst1), &
         wfir2(:,:,jst2),wfir1(:,:,jst1),zrhomt(:,:,:,j1),zrhoir(:,j1))
        zcc(ig,j1,j2)=zfinp(.false.,zrhomt(:,:,:,j1),expgqmt,zrhoir(:,j1), &
         expgir(:,ig))
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
! matrix elements for full BSE kernel if required
    if (bsefull) then
! compute the <v|exp(i(G+q).r)|c'> matrix elements
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,jst2,j2)
!$OMP DO
      do i1=1,nvbse
        ist1=istbse(i1,ik1)
        do j2=1,ncbse
          jst2=jstbse(j2,ik2)
          call genzrho(.false.,wfmt2(:,:,:,:,jst2),wfmt1(:,:,:,:,ist1), &
           wfir2(:,:,jst2),wfir1(:,:,ist1),zrhomt(:,:,:,i1),zrhoir(:,i1))
          zvc(ig,i1,j2)=zfinp(.false.,zrhomt(:,:,:,i1),expgqmt,zrhoir(:,i1), &
           expgir(:,ig))
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
! compute the <c|exp(i(G+q).r)|v'> matrix elements
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jst1,ist2,i2)
!$OMP DO
      do j1=1,ncbse
        jst1=jstbse(j1,ik1)
        do i2=1,nvbse
          ist2=istbse(i2,ik2)
          call genzrho(.false.,wfmt2(:,:,:,:,ist2),wfmt1(:,:,:,:,jst1), &
           wfir2(:,:,ist2),wfir1(:,:,jst1),zrhomt(:,:,:,j1),zrhoir(:,j1))
          zcv(ig,j1,i2)=zfinp(.false.,zrhomt(:,:,:,j1),expgqmt,zrhoir(:,j1), &
           expgir(:,ig))
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
    end if
! end loop over G-vectors
  end do
! get RPA inverse epsilon from file
  call getcf2pt(fname,vl,ngrpa,nwrpa,epsinv)
  t0=fourpi*wkptnr/omega
  do i1=1,nvbse
    do j1=1,ncbse
      a1=ijkbse(i1,j1,ik1)
      do i2=1,nvbse
        do j2=1,ncbse
          a2=ijkbse(i2,j2,ik2)
          zsum=0.d0
          do ig=1,ngrpa
            do jg=1,ngrpa
              t1=gqc(ig)*gqc(jg)
              if (t1.gt.1.d-6) then
                t1=t0/t1
                zsum=zsum+t1*epsinv(ig,jg,1)*zcc(ig,j1,j2)*conjg(zvv(jg,i1,i2))
              end if
            end do
          end do
          hmlbse(a1,a2)=hmlbse(a1,a2)-zsum
! compute off-diagonal blocks if required
          if (bsefull) then
            b1=a1+nbbse
            b2=a2+nbbse
            hmlbse(b1,b2)=hmlbse(b1,b2)+conjg(zsum)
            zsum=0.d0
            do ig=1,ngrpa
              do jg=1,ngrpa
                t1=gqc(ig)*gqc(jg)
                if (t1.gt.1.d-8) then
                  t1=t0/t1
                  zsum=zsum+t1*epsinv(ig,jg,1)*zcv(ig,j1,i2) &
                   *conjg(zvc(jg,i1,j2))
                end if
              end do
            end do
            hmlbse(a1,b2)=hmlbse(a1,b2)-zsum
            hmlbse(b1,a2)=hmlbse(b1,a2)+conjg(zsum)
          end if
! end loop over i2 and j2
        end do
      end do
! end loop over i1 and j1
    end do
  end do
! end loop over ik1
end do
deallocate(gqc,wfmt1,wfmt2,wfir1,wfir2,expqmt,expgqmt)
deallocate(zrhomt,zrhoir,zvv,zcc,epsinv)
if (bsefull) deallocate(zvc,zcv)
return
end subroutine
