
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmldbsek(ik2)
use modmain
implicit none
! arguments
integer, intent(in) :: ik2
! local variables
integer ik1,ig,jg
integer is,ia,ias,irc
integer i1,i2,j1,j2,a1,a2,b1,b2
integer ist1,ist2,jst1,jst2
real(8) vl(3),vc(3)
real(8) vgqc(3),t0,t1
complex(8) zsum
! allocatable arrays
real(8), allocatable :: gqc(:)
complex(8), allocatable :: wfmt1(:,:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:)
complex(8), allocatable :: wfir2(:,:,:)
complex(8), allocatable :: expqmt(:,:,:)
complex(8), allocatable :: expgqmt(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvv(:,:,:)
complex(8), allocatable :: zcc(:,:,:)
complex(8), allocatable :: zvc(:,:,:)
complex(8), allocatable :: zcv(:,:,:)
complex(8), allocatable :: epsinv(:,:,:)
! external functions
complex(8) zfinp
external zfinp
allocate(gqc(ngrpa))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(expqmt(lmmaxvr,nrcmtmax,natmtot))
allocate(expgqmt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvv(ngrpa,nvbse,nvbse))
allocate(zcc(ngrpa,ncbse,ncbse))
allocate(epsinv(nwrpa,ngrpa,ngrpa))
if (bsefull) then
  allocate(zvc(ngrpa,nvbse,ncbse))
  allocate(zcv(ngrpa,ncbse,nvbse))
end if
! generate the wavefunctions at k-point ik2
call genwfsvp(.false.,.false.,vkl(:,ik2),wfmt2,ngrtot,wfir2)
! begin loop over ik1
do ik1=1,nkptnr
! generate the wavefunctions at k-point ik1
  call genwfsvp(.false.,.false.,vkl(:,ik1),wfmt1,ngrtot,wfir1)
! q vector in lattice and Cartesian coordinates
  vl(:)=vkl(:,ik1)-vkl(:,ik2)
  vc(:)=vkc(:,ik1)-vkc(:,ik2)
! generate the function exp(iq.r) in the muffin-tins
  call genexpmt(vc,expqmt)
! loop over RPA G vectors
  do ig=1,ngrpa
! G+q vector in Cartesian coordinates
    vgqc(:)=vgc(:,ig)+vc(:)
! length of G+q vector
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
    do i1=1,nvbse
      ist1=istbse(i1,ik1)
      do i2=1,nvbse
        ist2=istbse(i2,ik2)
! note that the complex conjugate of the density is found because zfinp
! conjugates the first function
        call genzrho(.false.,wfmt2(:,:,:,:,ist2),wfmt1(:,:,:,:,ist1), &
         wfir2(:,:,ist2),wfir1(:,:,ist1),zrhomt,zrhoir)
        zvv(ig,i1,i2)=zfinp(.false.,zrhomt,expgqmt,zrhoir,expgir(:,ig))
      end do
    end do
! compute the <c|exp(i(G+q).r)|c'> matrix elements
    do j1=1,ncbse
      jst1=jstbse(j1,ik1)
      do j2=1,ncbse
        jst2=jstbse(j2,ik2)
        call genzrho(.false.,wfmt2(:,:,:,:,jst2),wfmt1(:,:,:,:,jst1), &
         wfir2(:,:,jst2),wfir1(:,:,jst1),zrhomt,zrhoir)
        zcc(ig,j1,j2)=zfinp(.false.,zrhomt,expgqmt,zrhoir,expgir(:,ig))
      end do
    end do
! matrix elements for full BSE kernel if required
    if (bsefull) then
! compute the <v|exp(i(G+q).r)|c'> matrix elements
      do i1=1,nvbse
        ist1=istbse(i1,ik1)
        do j2=1,ncbse
          jst2=jstbse(j2,ik2)
          call genzrho(.false.,wfmt2(:,:,:,:,jst2),wfmt1(:,:,:,:,ist1), &
           wfir2(:,:,jst2),wfir1(:,:,ist1),zrhomt,zrhoir)
          zvc(ig,i1,j2)=zfinp(.false.,zrhomt,expgqmt,zrhoir,expgir(:,ig))
        end do
      end do
! compute the <c|exp(i(G+q).r)|v'> matrix elements
      do j1=1,ncbse
        jst1=jstbse(j1,ik1)
        do i2=1,nvbse
          ist2=istbse(i2,ik2)
          call genzrho(.false.,wfmt2(:,:,:,:,ist2),wfmt1(:,:,:,:,jst1), &
           wfir2(:,:,ist2),wfir1(:,:,jst1),zrhomt,zrhoir)
          zcv(ig,j1,i2)=zfinp(.false.,zrhomt,expgqmt,zrhoir,expgir(:,ig))
        end do
      end do
    end if
! end loop over G vectors
  end do
! get RPA inverse epsilon from file
  call getepsinv_rpa(vl,epsinv)
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
                zsum=zsum+t1*epsinv(1,ig,jg)*zcc(ig,j1,j2)*conjg(zvv(jg,i1,i2))
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
                if (t1.gt.1.d-6) then
                  t1=t0/t1
                  zsum=zsum+t1*epsinv(1,ig,jg)*zcv(ig,j1,i2) &
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
