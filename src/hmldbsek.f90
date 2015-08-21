
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmldbsek(ik2)
use modmain
implicit none
! arguments
integer, intent(in) :: ik2
! local variables
integer ik1,ngknr,igk,ig,jg
integer is,ia,ias,irc
integer i1,i2,j1,j2,a1,a2
integer ist1,ist2,jst1,jst2
real(8) vl(3),vc(3)
real(8) vgqc(3),t0,t1
complex(8) zsum
! allocatable arrays
integer, allocatable :: igkignr(:)
real(8), allocatable :: vgklnr(:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: gkcnr(:)
real(8), allocatable :: tpgkcnr(:,:)
real(8), allocatable :: gqc2(:)
complex(8), allocatable :: sfacgknr(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
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
complex(8), allocatable :: epsinv(:,:,:)
! external functions
complex(8) zfinp
external zfinp
allocate(igkignr(ngkmax))
allocate(vgklnr(3,ngkmax))
allocate(vgkcnr(3,ngkmax))
allocate(gkcnr(ngkmax))
allocate(tpgkcnr(2,ngkmax))
allocate(gqc2(ngrpa))
allocate(sfacgknr(ngkmax,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
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
! generate G+k vectors for k-point ik2
call gengpvec(vkl(:,ik2),vkc(:,ik2),ngknr,igkignr,vgklnr,vgkcnr)
! generate the spherical coordinates of the G+k vectors
do igk=1,ngknr
  call sphcrd(vgkcnr(:,igk),gkcnr(igk),tpgkcnr(:,igk))
end do
! generate the structure factors
call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
! get the eigenvector for k-point ik2
call getevecfv(vkl(:,ik2),vgklnr,evecfv)
call getevecsv(vkl(:,ik2),evecsv)
! generate the wavefunctions for k-point ik2
call genwfsv(.false.,.false.,.false.,ngknr,igkignr,evalsv,apwalm,evecfv, &
 evecsv,wfmt2,ngrtot,wfir2)
! begin loop over ik1
do ik1=1,nkptnr
! generate G+k vectors for k-point ik1
  call gengpvec(vkl(:,ik1),vkc(:,ik1),ngknr,igkignr,vgklnr,vgkcnr)
! generate the spherical coordinates of the G+k vectors
  do igk=1,ngknr
    call sphcrd(vgkcnr(:,igk),gkcnr(igk),tpgkcnr(:,igk))
  end do
! generate the structure factors
  call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
  call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
! get the eigenvectors from file
  call getevecfv(vkl(:,ik1),vgklnr,evecfv)
  call getevecsv(vkl(:,ik1),evecsv)
! generate the wavefunctions for k-point ik1
  call genwfsv(.false.,.false.,.false.,ngknr,igkignr,evalsv,apwalm,evecfv, &
   evecsv,wfmt1,ngrtot,wfir1)
! q vector in lattice and Cartesian coordinates
  vl(:)=vkl(:,ik1)-vkl(:,ik2)
  vc(:)=vkc(:,ik1)-vkc(:,ik2)
! generate the function exp(iq.r) in the muffin-tins
  call genexpmt(vc,expqmt)
! loop over G vectors
  do ig=1,ngrpa
! G+q vector in Cartesian coordinates
    vgqc(:)=vgc(:,ig)+vc(:)
! length of G+q vector squared
    gqc2(ig)=vgqc(1)**2+vgqc(2)**2+vgqc(3)**2
! compute the fuction exp(i(G+q).r) in the muffin-tins
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
          do jg=1,ngrpa
            t1=t0/(gqc2(jg)+lmda2rpa)
            do ig=1,ngrpa
              zsum=zsum+t1*epsinv(1,ig,jg)*zcc(ig,j1,j2)*conjg(zvv(jg,i1,i2))
            end do
          end do
          hmlbse(a1,a2)=hmlbse(a1,a2)-zsum
        end do
      end do
    end do
  end do
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr,gqc2)
deallocate(sfacgknr,apwalm,evecfv,evecsv)
deallocate(wfmt1,wfmt2,wfir1,wfir2,expqmt,expgqmt)
deallocate(zrhomt,zrhoir,zvv,zcc,epsinv)
return
end subroutine
