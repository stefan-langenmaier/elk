
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phonon
use modmain
use modphonon
implicit none
integer nsym,isym,sym(3,3,48)
integer nppt,npptnr,ik,ip,ispn
integer igp
real(8) wpptnr
real(8) s(3,3),vl(3),vc(3),t1
! allocatable arrays
integer, allocatable :: ipmap(:,:,:)
integer, allocatable :: ipmapnr(:,:,:)
integer, allocatable :: ivp(:,:)
integer, allocatable :: ngp(:,:)
integer, allocatable :: ngpq(:,:)
integer, allocatable :: igpig(:,:,:)
integer, allocatable :: igpqig(:,:,:)
real(8), allocatable :: vpl(:,:)
real(8), allocatable :: vpc(:,:)
real(8), allocatable :: wppt(:)
real(8), allocatable :: vgpl(:,:,:,:)
real(8), allocatable :: vgpql(:,:,:,:)
real(8), allocatable :: vgpc(:,:,:,:)
real(8), allocatable :: vgpqc(:,:,:,:)
real(8), allocatable :: gpc(:,:,:)
real(8), allocatable :: gpqc(:,:,:)
real(8), allocatable :: tpgpc(:,:,:,:)
real(8), allocatable :: tpgpqc(:,:,:,:)
complex(8), allocatable :: sfacgp(:,:,:,:)
complex(8), allocatable :: sfacgpq(:,:,:,:)
complex(8), allocatable :: dwfpw(:,:,:)
! initialise universal variables
call init0
call init1
call init2
! allocate local arrays
allocate(ipmap(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
npptnr=ngridk(1)*ngridk(2)*ngridk(3)
allocate(ivp(3,npptnr))
allocate(ngp(nspnfv,npptnr))
allocate(ngpq(nspnfv,npptnr))
allocate(igpig(ngkmax,nspnfv,npptnr))
allocate(igpqig(ngkmax,nspnfv,npptnr))
allocate(vpl(3,npptnr))
allocate(vpc(3,npptnr))
allocate(wppt(npptnr))
allocate(vgpl(3,ngkmax,nspnfv,npptnr))
allocate(vgpql(3,ngkmax,nspnfv,npptnr))
allocate(vgpc(3,ngkmax,nspnfv,npptnr))
allocate(vgpqc(3,ngkmax,nspnfv,npptnr))
allocate(gpc(ngkmax,nspnfv,npptnr))
allocate(gpqc(ngkmax,nspnfv,npptnr))
allocate(tpgpc(2,ngkmax,nspnfv,npptnr))
allocate(tpgpqc(2,ngkmax,nspnfv,npptnr))
allocate(sfacgp(ngkmax,natmtot,nspnfv,npptnr))
allocate(sfacgpq(ngkmax,natmtot,nspnfv,npptnr))
allocate(dwfpw(ngkmax,nspinor,nstsv))
!***** deallocate
! read density and potential from file
call readstate
! read in the eigenvalues and occupancies
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
10 continue
call dyntask(80)
if (iqph.eq.0) return
! phonon dry run
if (task.eq.201) then
  close(80)
  goto 10
end if
! find all symmetries which leave q invariant (little group)
nsym=0
do isym=1,nsymkpt
  s(:,:)=dble(symkpt(:,:,isym))
  call r3mtv(s,vql(:,iqph),vl)
  t1=abs(vql(1,iqph)-vl(1))+abs(vql(2,iqph)-vl(2))+abs(vql(3,iqph)-vl(3))
  if (t1.lt.epslat) then
    nsym=nsym+1
    sym(:,:,nsym)=symkpt(:,:,isym)
  end if
end do
! generate the reduced k-point set for current q vector perturbation
call genppts(.false.,nsym,sym,ngridk,npptnr,epslat,bvec,kptboxl,nppt,ipmap, &
 ipmapnr,ivp,vpl,vpc,wppt,wpptnr)
do ip=1,nppt
  do ispn=1,nspnfv
    vl(:)=vpl(:,ip)
    vc(:)=vpc(:,ip)
! spin-spiral case
    if (spinsprl) then
      if (ispn.eq.1) then
        vl(:)=vl(:)+0.5d0*vqlss(:)
        vc(:)=vc(:)+0.5d0*vqcss(:)
      else
        vl(:)=vl(:)-0.5d0*vqlss(:)
        vc(:)=vc(:)-0.5d0*vqcss(:)
      end if
    end if
! generate the G+k vectors
    call gengpvec(vl,vc,ngp(ispn,ip),igpig(:,ispn,ip),vgpl(:,:,ispn,ip), &
     vgpc(:,:,ispn,ip))
! generate the spherical coordinates of the G+k vectors
    do igp=1,ngp(ispn,ip)
      call sphcrd(vgpc(:,igp,ispn,ip),gpc(igp,ispn,ip),tpgpc(:,igp,ispn,ip))
    end do
! generate structure factors for G+k vectors
    call gensfacgp(ngp(ispn,ip),vgpc(:,:,ispn,ip),ngkmax,sfacgp(:,:,ispn,ip))
! k+q vector
    vl(:)=vl(:)+vql(:,iqph)
    vc(:)=vc(:)+vqc(:,iqph)
! generate the G+k+q vectors, etc.
    call gengpvec(vl,vc,ngpq(ispn,ip),igpqig(:,ispn,ip),vgpql(:,:,ispn,ip), &
     vgpqc(:,:,ispn,ip))
    do igp=1,ngpq(ispn,ip)
      call sphcrd(vgpqc(:,igp,ispn,ip),gpqc(igp,ispn,ip),tpgpqc(:,igp,ispn,ip))
    end do
    call gensfacgp(ngpq(ispn,ip),vgpqc(:,:,ispn,ip),ngkmax,sfacgpq(:,:,ispn,ip))
  end do
end do
! begin the self-consistent loop
do iscl=1,maxscl
! zero the density and magnetisation derivatives
  drhomt(:,:,:)=0.d0
  drhoir(:)=0.d0
  if (spinpol) then
    dmagmt(:,:,:,:)=0.d0
    dmagir(:,:)=0.d0
  end if
! loop over k-points
  do ip=1,nppt
! compute the first-order change in the wavefunction
    call phdwfpw(vpl(:,ip),dwfpw)
! add to the density and magnetisation derivatives
    call drhomagk(vpl(:,ip),ngp(:,ip),ngpq(:,ip),igpig(:,:,ip),igpqig(:,:,ip), &
     vgpl(:,:,:,ip),gpqc(:,:,ip),tpgpqc(:,:,:,ip),sfacgpq(:,:,:,ip),dwfpw)
  end do
! convert muffin-tin density/magnetisation derivatives to spherical harmonics
  call drhomagsh
! end the self-consistent loop
end do
goto 10
return
end subroutine

