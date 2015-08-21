
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ephcouple
use modmain
use modphonon
use modmpi
use modstore
implicit none
! local variables
integer iq,ik,jk,jkq
integer iv(3),ist,jst
integer is,ia,ias,js,ja,jas
integer ip,ir,irc,isym,i,j,n
real(8) vl(3),x,t1,t2,t3,t4,t5
complex(8) z1
! allocatable arrays
real(8), allocatable :: wq(:,:),gq(:,:)
complex(8), allocatable :: dynq(:,:,:),ev(:,:)
complex(8), allocatable :: dvphmt(:,:,:,:),dvphir(:,:)
complex(8), allocatable :: zfmt(:,:),gzfmt(:,:,:,:)
complex(8), allocatable :: ephmat(:,:,:)
! external functions
real(8) sdelta,stheta
external sdelta,stheta
! increase the angular momentum cut-off on the inner part of the muffin-tin
lmaxinr0=lmaxinr
lmaxinr=max(lmaxinr,4)
! initialise universal variables
call init0
call init1
call init2
! check k-point grid is commensurate with q-point grid
iv(:)=mod(ngridk(:),ngridq(:))
if ((iv(1).ne.0).or.(iv(2).ne.0).or.(iv(3).ne.0)) then
  write(*,*)
  write(*,'("Error(ephcouple): k-point grid incommensurate with q-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" ngridq : ",3I6)') ngridq
  write(*,*)
  stop
end if
! allocate global arrays
if (allocated(dvsmt)) deallocate(dvsmt)
allocate(dvsmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(dvsir)) deallocate(dvsir)
allocate(dvsir(ngtot))
! allocate local arrays
allocate(wq(nbph,nqpt),gq(nbph,nqpt))
allocate(dynq(nbph,nbph,nqpt),ev(nbph,nbph))
allocate(dvphmt(lmmaxvr,nrcmtmax,natmtot,nbph))
allocate(dvphir(ngtot,nbph))
allocate(zfmt(lmmaxvr,nrmtmax))
allocate(gzfmt(lmmaxvr,nrmtmax,3,natmtot))
! read in the density and potentials from file
call readstate
! read in the Fermi energy
call readfermi
! find the linearisation energies
call linengy
! set the speed of light >> 1 (non-relativistic approximation)
solsc=sol*100.d0
! new file extension for eigenvector files with c >> 1
filext='_EPH.OUT'
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! generate the spin-orbit coupling radial functions
call gensocfr
! generate the first- and second-variational eigenvectors and eigenvalues
call genevfsv
! restore the speed of light
solsc=sol
! compute the occupancies and density of states at the Fermi energy
call occupy
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! loop over all atoms
do ias=1,natmtot
  is=idxis(ias)
! convert potential to complex spherical harmonic expansion
  call rtozfmt(nrmt(is),nrmtinr(is),1,vsmt(:,:,ias),1,zfmt)
! compute the gradients of the Kohn-Sham potential for the rigid-ion term
  call gradzfmt(nrmt(is),nrmtinr(is),spr(:,is),zfmt,nrmtmax,gzfmt(:,:,:,ias))
end do
! loop over phonon q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(ephcouple): ",I6," of ",I6," q-points")') iq,nqpt
! diagonalise the dynamical matrix
  call dynev(dynq(:,:,iq),wq(:,iq),ev)
! loop over phonon branches
  do j=1,nbph
! zero any negative frequencies
    if (wq(j,iq).lt.0.d0) wq(j,iq)=0.d0
! find change in Kohn-Sham potential for mode j
    dvphmt(:,:,:,j)=0.d0
    dvphir(:,j)=0.d0
    i=0
    do is=1,nspecies
! prefactor
      t1=2.d0*spmass(is)*wq(j,iq)
      if (t1.gt.1.d-8) then
        t1=1.d0/sqrt(t1)
      else
        t1=0.d0
      end if
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ip=1,3
          i=i+1
! read in the Cartesian change in Kohn-Sham potential
          call readdvs(iq,is,ia,ip)
! add the rigid-ion term
          do ir=1,nrmt(is)
            dvsmt(:,ir,ias)=dvsmt(:,ir,ias)-gzfmt(:,ir,ip,ias)
          end do
! multiply with eigenvector component and add to total phonon potential
          z1=t1*ev(i,j)
          do js=1,nspecies
            do ja=1,natoms(js)
              jas=idxas(ja,js)
              irc=0
              do ir=1,nrmt(js),lradstp
                irc=irc+1
                dvphmt(:,irc,jas,j)=dvphmt(:,irc,jas,j)+z1*dvsmt(:,ir,jas)
              end do
            end do
          end do
          dvphir(:,j)=dvphir(:,j)+z1*dvsir(:)
        end do
      end do
    end do
  end do
! zero the phonon linewidths array
  gq(:,iq)=0.d0
! begin parallel loop over non-reduced k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ephmat,jk,vl,isym,jkq) &
!$OMP PRIVATE(t1,t2,t3,t4,t5,ist,jst,x,i)
!$OMP DO
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
    allocate(ephmat(nstsv,nstsv,nbph))
! equivalent reduced k-point
    jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! compute the electron-phonon coupling matrix elements
    call genephmat(iq,ik,dvphmt,dvphir,ephmat)
! k+q-vector in lattice coordinates
    vl(:)=vkl(:,ik)+vql(:,iq)
! index to k+q-vector
    call findkpt(vl,isym,jkq)
    t1=twopi*wkptnr*occmax
    do ist=1,nstsv
      x=(efermi-evalsv(ist,jkq))/swidth
      t2=stheta(stype,x)
      do jst=1,nstsv
        x=(efermi-evalsv(jst,jk))/swidth
        t3=stheta(stype,x)
! loop over phonon branches
        do i=1,nbph
          x=(wq(i,iq)+evalsv(jst,jk)-evalsv(ist,jkq))/swidth
          t4=sdelta(stype,x)/swidth
          t5=dble(ephmat(ist,jst,i))**2+aimag(ephmat(ist,jst,i))**2
!$OMP CRITICAL
          gq(i,iq)=gq(i,iq)+t1*(t3-t2)*t4*t5
!$OMP END CRITICAL
        end do
      end do
    end do
    deallocate(ephmat)
! end loop over k-points
  end do
!$OMP END DO
!$OMP END PARALLEL
! end loop over phonon q-points
end do
! add gq from each process and redistribute
if (np_mpi.gt.1) then
  n=nbph*nqpt
  call mpi_allreduce(mpi_in_place,gq,n,mpi_double_precision,mpi_sum, &
   mpi_comm_kpt,ierror)
end if
filext='.OUT'
if (mp_mpi) then
! write the phonon linewidths to file
  call writegamma(gq)
! write electron-phonon coupling constants to file
  call writelambda(wq,gq)
end if
deallocate(wq,gq,dynq,ev)
deallocate(dvphmt,dvphir)
deallocate(zfmt,gzfmt)
! restore lmaxinr
lmaxinr=lmaxinr0
return
end subroutine

