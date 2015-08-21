
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dielectricq
use modmain
implicit none
! local variables
integer iq,ik,jk,jkq
integer ngknr,ngkq
integer igk,igkq,isym
integer ist,jst,iw
real(8) vkql(3),vkqc(3)
real(8) v(3),eji,q2,t1
complex(8) zsum,zt1,zt2
character(256) fname
! allocatable arrays
integer, allocatable :: igkignr(:)
integer, allocatable :: igkqig(:)
integer, allocatable :: map(:)
real(8), allocatable :: vgklnr(:,:)
real(8), allocatable :: vgkql(:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: vgkqc(:,:)
real(8), allocatable :: gkcnr(:)
real(8), allocatable :: gkqc(:)
real(8), allocatable :: tpgkcnr(:,:)
real(8), allocatable :: tpgkqc(:,:)
real(8), allocatable :: w(:)
complex(8), allocatable :: wfpwk(:,:,:)
complex(8), allocatable :: wfpwkq(:,:,:)
complex(8), allocatable :: chi(:)
complex(8), allocatable :: eps(:)
! initialise universal variables
call init0
call init1
call init2
! allocate local arrays
allocate(igkignr(ngkmax))
allocate(igkqig(ngkmax))
allocate(map(ngkmax))
allocate(vgklnr(3,ngkmax))
allocate(vgkql(3,ngkmax))
allocate(vgkcnr(3,ngkmax))
allocate(vgkqc(3,ngkmax))
allocate(gkcnr(ngkmax))
allocate(gkqc(ngkmax))
allocate(tpgkcnr(2,ngkmax))
allocate(tpgkqc(2,ngkmax))
allocate(w(nwdos))
allocate(wfpwk(ngkmax,nspinor,nstsv))
allocate(wfpwkq(ngkmax,nspinor,nstsv))
allocate(chi(nwdos))
allocate(eps(nwdos))
! output the q-point set to file
call writeqpts
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! read in the eigenvalues and occupancies
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! generate energy grid (starting from zero)
t1=wdos(2)/dble(nwdos)
do iw=1,nwdos
  w(iw)=t1*dble(iw-1)
end do
do iq=1,nqpt
  q2=vqc(1,iq)**2+vqc(2,iq)**2+vqc(3,iq)**2
  if (sqrt(q2).lt.epslat) goto 20
  write(*,'("Info(dielectricq): ",I6," of ",I6," q-points")') iq,nqpt
  chi(:)=0.d0
  do ik=1,nkptnr
! equivalent reduced k-point
    jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
! generate G+k vectors
    call gengpvec(vklnr(:,ik),vkcnr(:,ik),ngknr,igkignr,vgklnr,vgkcnr,gkcnr, &
     tpgkcnr)
! get the plane wave wavefunction at k
    call getwfpw(vklnr(:,ik),vgklnr,wfpwk)
! k+q vector in lattice coordinates
    vkql(:)=vklnr(:,ik)+vql(:,iq)
! find equivalent reduced k-point
    call findkpt(vkql,isym,jkq)
! k+q in Cartesian coordinates
    call r3mv(bvec,vkql,vkqc)
! generate G+k+q vectors
    call gengpvec(vkql,vkqc,ngkq,igkqig,vgkql,vgkqc,gkqc,tpgkqc)
! get the plane wave wavefunction at k+q
    call getwfpw(vkql,vgkql,wfpwkq)
! generate the map from G+k to G+k+q
    map(:)=0
    do igk=1,ngknr
      v(:)=vgklnr(:,igk)+vql(:,iq)
      do igkq=1,ngkq
        t1=abs(v(1)-vgkql(1,igkq)) &
          +abs(v(2)-vgkql(2,igkq)) &
          +abs(v(3)-vgkql(3,igkq))
        if (t1.lt.epslat) then
          map(igk)=igkq
          goto 10
        end if
      end do
10 continue
    end do
    do ist=1,nstsv
      do jst=1,nstsv
        eji=evalsv(jst,jkq)-evalsv(ist,jk)
        zsum=0.d0
        do igk=1,ngknr
          igkq=map(igk)
          if (igkq.ne.0) then
            zsum=zsum+conjg(wfpwkq(igkq,1,ist))*wfpwk(igk,1,jst)
            if (spinpol) zsum=zsum+conjg(wfpwkq(igkq,2,ist))*wfpwk(igk,2,jst)
          end if
        end do
        t1=dble(zsum)**2+aimag(zsum)**2
        t1=t1*occsv(ist,jk)*(1.d0-occsv(jst,jkq)/occmax)
        t1=t1*wkptnr(ik)
        do iw=1,nwdos
          zt1=w(iw)-eji+zi*swidth
          zt2=w(iw)+eji-zi*swidth
          chi(iw)=chi(iw)+t1*(1.d0/zt1-1.d0/zt2)
        end do
      end do
    end do
  end do
  chi(:)=(fourpi/omega)*chi(:)
! inverse dielectric function
  eps(:)=1.d0+chi(:)/q2
  write(fname,'("EPSILON_Q",I4.4,".OUT")') iq
  open(60,file=trim(fname),action='WRITE',form='FORMATTED')
  do iw=1,nwdos
    write(60,'(2G18.10)') w(iw),dble(eps(iw))
  end do
  write(60,'("     ")')
  do iw=1,nwdos
    write(60,'(2G18.10)') w(iw),aimag(eps(iw))
  end do
  close(60)
20 continue
end do
deallocate(igkignr,igkqig,map)
deallocate(vgklnr,vgkql,vgkcnr,vgkqc)
deallocate(gkcnr,gkqc,tpgkcnr,tpgkqc,w)
deallocate(wfpwk,wfpwkq,chi,eps)
return
end subroutine

