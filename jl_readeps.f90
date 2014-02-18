!===============================================================================
!
! Utilities:
!
! (1) JLreadeps    
!
! Serial code that reads epsmat
!
! epsinvomega.inp:
! eps0mat      ! epsmat file, full-frequency or static
!===============================================================================

#include "f_defs.h"

program jl_readeps

  use global_m
  implicit none

  integer :: i,j,k,iq,mq,ig,igp,itape,indx,jj,np,nqnonzero
  real(DP) :: omega
  
  character*256, parameter :: fninp = "JLepsinvomega.inp"
  character*256 :: fneps,fnwfn,fnrho,fngpp,fnffr,fnffa
  real(DP) :: q(3)
  integer :: g(3),gp(3),gmgp(3),nband,ktot,ntranq
  
  integer :: freq_dep,nFreq,ii,nq,ng,nmtx,kgrid(3),kmax(3)
  real(DP) :: dDeltaFreq,dBrdning,ecuts,delta,qvec(3)
  real(DP), allocatable :: qpt(:,:)
  real(DP), allocatable :: dFreqGrid(:),ekin(:)
  integer, allocatable :: gvec(:,:)
  integer, allocatable :: isrtx(:)
  integer, allocatable :: isorti(:)
  SCALAR, allocatable :: eps(:)
  complex(DPC), allocatable :: dFreqBrd(:)
  complex(DPC), allocatable :: epsPP(:)
  complex(DPC), allocatable :: epsR(:)
  complex(DPC), allocatable :: epsA(:)
  
  real(DP) :: bdot(3,3),celvol,ebind
  integer :: nvecs,nspin
  
  integer nproc_para,num_gvec,gx,gy,gz
  complex(DPC) :: epsStatic,xcdum(2)
  real(DP) :: rho0
  SCALAR :: rhogmgp
  
  real(DP) :: wp2,qg(3),qgp(3),qgqg,qgqgp,lambda,phi
  SCALAR :: Omega2,wtilde2,epsggp,I_epsggp,eps_static,eps_dynamic
  complex(DPC) :: wtilde2_temp

!-----------------------------
  call open_file(unit=9,file='fulleps.dat',form='formatted',status='replace')
  call open_file(unit=8,file='gvecs.dat',form='formatted',status='replace')

! read input file
  write(6,'(/,1x,"reading",1x,a,1x,"file",/)')trim(fninp)  
  call open_file(55,file=trim(fninp),form='formatted',status='old')
  read(55,'(a)') fneps  
  call close_file(55)
  
!-----------------------------
! read eps file
  write(6,'(1x,"reading",1x,a,1x,"file",/)')trim(fneps)

  itape=12
  call open_file(unit=itape,file=trim(fneps),form='unformatted',status='old')

  read(itape)
  read(itape) freq_dep,ii
  write(6,*) 'Frequency dependence ', freq_dep
  if (freq_dep.eq.2) then
     nFreq=ii
     write(6,*) 'Number of Fequencies ', nFreq
     SAFE_ALLOCATE(epsR, (nFreq))
     SAFE_ALLOCATE(dFreqGrid,(nFreq))
     SAFE_ALLOCATE(dFreqBrd,(nFreq))
  endif
  
  read(itape) (kgrid(i),i=1,3)

  if (freq_dep.eq.2) then
     read(itape) (dFreqGrid(i),i=1,nFreq),(dFreqBrd(i),i=1,nFreq)
     if (nFreq.gt.1) dDeltaFreq=dFreqGrid(2)-dFreqGrid(1)
     dBrdning=IMAG(dFreqBrd(1))
     write(6,*) 'broadening=', dBrdning
  else
     read(itape)
  endif
  read(itape)
  read(itape)
  read(itape) ecuts
  write(6,*) 'Screened Coulomb cutoff', ecuts

  read(itape) nqnonzero
  write(6,*) 'nqnonzero', nqnonzero
  nq = nqnonzero

  read(itape) ng
  write(6,*) 'Number of Gvecs', ng
  !read(itape) ng,ktot,kmax(1:3)
  !write(6,*) 'Number of G-vectors, Number of k points', ng,ktot
  !backspace(itape)

  !SAFE_ALLOCATE(gvec,(3,ng))
  !SAFE_ALLOCATE(ekin,(ng))
  !SAFE_ALLOCATE(isrtx,(ng))

  !read(itape) ng,ktot,kmax(1:3),((gvec(jj,ig),jj=1,3),ig=1,ng)

  !read(itape) nq
  !write(6,*) 'Number of q vectors', nq
  !if( trim(fneps) .eq. "chi0mat") nq = 1

  do ii = 1,nq

     read(itape) ng,nmtx! epsinv.f90
     write(6,*) 'nmtx', nmtx
     read(itape) ! epsinv.f90
     read(itape) ! epsinv.f90


     !read(itape) ntranq
     !write(6,*) 'ntranq ', ntranq

     !read(itape) nmtx,np,(isrtx(jj),ekin(jj),jj=1,ng)
     !write(6,*) 'Size of epsilon ', nmtx

     !do jj = 1,nmtx
     !   write(8,*)gvec(1:3,isrtx(jj))
     !enddo
     
     do j=1,nmtx
        do i=1,nmtx
           read(itape) (epsR(k),k=1,nFreq)
           do indx = 1,nFreq
              write(9,100)epsR(indx)
           enddo
        enddo
#ifdef CPLX
        do i=1,nmtx
           read(itape)
        enddo
#endif
     enddo
  enddo

  
  call close_file(itape)
  
  write(6,'(a,i6)')     "     omega num  = ", nFreq
  call close_file(9)
  call close_file(8)

!-----------------------------
! deallocate and finish
  
  SAFE_DEALLOCATE(epsPP)
  SAFE_DEALLOCATE(epsR)
  SAFE_DEALLOCATE(epsA)
  SAFE_DEALLOCATE(dFreqGrid)
  SAFE_DEALLOCATE(dFreqBrd)
  
100 format(3f25.15)
  
end program jl_readeps
