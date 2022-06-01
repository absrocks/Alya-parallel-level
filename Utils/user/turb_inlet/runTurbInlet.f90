program runTurbInlet

  ! INPUT:
  !   (1) inletTurb.dat // turbulent database generated with genTurbInlet.py 
  !
  ! OUTPUT:
  !   (1) pipe.ensi.VELAV-000001      // global index for inlet nodes => for Alya
  !   (2) pipe.ensi.VELFL-000001      // coordinates of index nodes   => for script
  !
  ! D. Mira, Last update: 19/12/2016

  implicit none

  integer  :: nstep,nboun,ndime
  integer  :: ii,jj,kk,iboun,istep
  real(kind=8)     :: auxi
  ! beware with real 4 there can be some errores
  real(kind=8), allocatable :: veloc(:,:,:),avvel(:,:),avve2(:,:)
  character(len=50) :: arg

  open(10,file='inletTurb.dat',status='unknown')
  open(30,file='inlet.ensi.VELAV-000001',status='unknown')
  open(40,file='inlet.ensi.VELFL-000001',status='unknown')

  nstep    =     401 ! # Temporal steps in the turbulent database  
!  nboun    =     15553! # Nodes of the turbulent database
  call getarg(1,arg)
  read(arg,*)nboun
  print*,'nboun',nboun
  ndime    =        3 ! # Velocity components
  !
  ! Allocation
  !
  allocate( veloc(nstep,nboun,ndime) )
  allocate( avvel(nboun,ndime) )
  allocate( avve2(nboun,ndime) )
  !
  ! Read in the turbulent database
  !
  veloc = 0.0
  do istep = 1,nstep
     do iboun = 1,nboun
        read(10,*)ii,jj,veloc(istep,iboun,1:ndime)
     end do
  end do
  !
  ! Time-averaging the turbulent database 
  !
  avvel = 0.0
  avve2 = 0.0
  do istep = 1,nstep
     do iboun = 1,nboun
        avvel(iboun,1:ndime) = avvel(iboun,1:ndime) + veloc(istep,iboun,1:ndime)
        avve2(iboun,1:ndime) = avve2(iboun,1:ndime) + veloc(istep,iboun,1:ndime)**2
     end do
!     write(778,'(a,i9,9(1x,e14.7))')'istep,veloc,avvel,avve2',istep,veloc(istep,1,:),avvel(1,:),avve2(1,:)
  end do

  avvel(1:nboun,1:ndime) = avvel(1:nboun,1:ndime) / (nstep*1.0)
  avve2(1:nboun,1:ndime) = avve2(1:nboun,1:ndime) / (nstep*1.0)

  !
  ! Generate <u_i> ensight file
  !
  write(30,'(A51)')'Alya Ensight Gold --- Vector per-node variables file'
  write(30,'(A4)' )'part'
  write(30,'(A10)')'         1'
  write(30,'(A11)')'coordinates'

  do ii=1,ndime
     do jj=1,nstep
        do kk=1,nboun
           write(30,*)avvel(kk,ii)
        end do
     end do
  end do
  !
  ! Generate <u'_i*u'_i> ensight file
  !
  write(40,'(A51)')'Alya Ensight Gold --- Vector per-node variables file'
  write(40,'(A4)' )'part'
  write(40,'(A10)')'         1'
  write(40,'(A11)')'coordinates'

  do ii=1,ndime
     do jj=1,nstep
        do kk=1,nboun
           auxi= avve2(kk,ii) - avvel(kk,ii)**2
           if (auxi<0.0) write(777,*)kk,avve2(kk,ii),avvel(kk,ii),auxi
           write(40,*) sqrt(max(auxi,0.0))
        end do
     end do
  end do

  deallocate(veloc)
  deallocate(avvel)
  deallocate(avve2)

  close(10)
  close(30)
  close(40)
  
  
end program

  ! --  Alya ensight format  --
  !
  !     u (ipoin)
  !     ...
  !     u (npoin)
  !     v (ipoin)
  !     ...
  !     v (npoin)
  !     w (ipoin)
  !     ...
  !     w (npoin)
  !


!!  do ii =1,ntime
!!     !ipass = int( mod(istep*dt,dt_turb*nstep)/dt_turb )
!!     !alpha = mod( istep*dt,dt_turb) / dt_turb
!!     do istep = 1,nstep
!!        do iboun = 1,nboun
!!           avvel(istep,iboun,1:ndime) = avvel(istep,iboun,1:ndime) + veloc(istep,iboun,1:ndime) 
!!        end do
!!     end do 
!!  end do
!!
!!  avvel(1:nstep,1:nboun,1:ndime) = avvel(1:nstep,1:nboun,1:ndime) / ntime


