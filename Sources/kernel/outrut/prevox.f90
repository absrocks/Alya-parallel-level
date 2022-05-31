subroutine prevox()
  !-----------------------------------------------------------------------
  !****
  ! NAME
  !    prevox
  ! DESCRIPTION
  !    This routine preprocesses information needed for saving 
  !    to voxel format.
  !    Each slave scans its own domain to see if any voxels lie inside it. It keeps a 
  !    registry of which voxels are inside (by index) and, more importantly, inside
  !    which element each voxel is. This is kept and used later to send voxel data.
  !    The second part of the routine is the slaves transmitting their voxel list
  !    to the master, so that afterwords it only needs to receive the voxel data 
  !    without info about its location.
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_postpr
  use mod_elsest,   only : elsest_host_element
  use mod_parall,   only : par_bin_comin
  use mod_parall,   only : par_bin_comax
  use mod_messages, only : livinf
  implicit none

  integer(ip)          :: ix,iy,iz,ielem,is
  integer(ip)          :: MyVoxelBox(2,3),VoxelGuess
  integer(ip), target  :: dummp(1)
  real(rp)             :: shapp(mnode),deriv(mnode*3),coloc(3)
  integer(ip),pointer  :: TempIndexList(:),TempElemList(:)
  integer(ip)          :: OldSize,kx,ky,kz
  real(rp)             :: Position(3),xx,yy,zz,rx,ry,rz,dista
 
  if( kfl_abovx == 0 ) return

  call livinf(0_ip,'CREATE VOXELS',0_ip)

  dummp = 0
  if( kfl_abovx == 1 ) then
     bobvx(1,1:ndime) = par_bin_comin(1:ndime)
     bobvx(2,1:ndime) = par_bin_comax(1:ndime)
  end if

  if (INOTSLAVE) then ! Master allocates memory for indexes of displacements and number of voxels
     allocate (nslvx(npart+1),dplvx(npart+1)) 
  endif

  if( INOTMASTER ) then
     !
     ! This is the spacing between voxel points
     !
     kx = max( resvx(1) , 2_ip )
     ky = max( resvx(2) , 2_ip )
     kz = max( resvx(3) , 2_ip )
     xx = ( bobvx(2,1)-bobvx(1,1) ) / real((kx-1),rp)  ! The -1 is for positioning from 0 to N-1
     yy = ( bobvx(2,2)-bobvx(1,2) ) / real((ky-1),rp)
     zz = ( bobvx(2,3)-bobvx(1,3) ) / real((kz-1),rp)
     !
     ! To speed up, check only the points that correspond to the inside of the domain
     !
     MyVoxelBox(1,1) = floor   ( ( xmima(1,1)- bobvx(1,1) ) / xx )
     MyVoxelBox(1,2) = floor   ( ( xmima(1,2)- bobvx(1,2) ) / yy )
     MyVoxelBox(1,3) = floor   ( ( xmima(1,3)- bobvx(1,3) ) / zz )
     MyVoxelBox(2,1) = ceiling ( ( xmima(2,1)- bobvx(1,1) ) / xx ) + 1
     MyVoxelBox(2,2) = ceiling ( ( xmima(2,2)- bobvx(1,2) ) / yy ) + 1
     MyVoxelBox(2,3) = ceiling ( ( xmima(2,3)- bobvx(1,3) ) / zz ) + 1
     VoxelGuess      =    ( MyVoxelBox(2,1) - MyVoxelBox(1,1) ) &
          &             * ( MyVoxelBox(2,2) - MyVoxelBox(1,2) ) &
          &             * ( MyVoxelBox(2,3) - MyVoxelBox(1,3) ) 
     if (VoxelGuess < 1) VoxelGuess = 1 ! Array must be at least size 1
     if (ISLAVE) then ! guess number of voxels inside domain
        allocate(TempIndexList(VoxelGuess),TempElemList(VoxelGuess))
     else ! We are sequential, allocate for all voxels just in case
        allocate(TempIndexList(product(resvx)),TempElemList(product(resvx)))
     endif
     !
     !Search voxels contained inside the domain
     !
     nvoxl   = 0
     OldSize = size(TempElemList)

     do iz= max(1_ip,MyVoxelBox(1,3)), min(resvx(3),MyVoxelBox(2,3))
        rz = real(iz-1,rp)
        do iy = max(1_ip,MyVoxelBox(1,2)), min(resvx(2),MyVoxelBox(2,2))
           ry = real(iy-1,rp)
           do  ix =  max(1_ip,MyVoxelBox(1,1)), min(resvx(1),MyVoxelBox(2,1))
              rx = real(ix-1,rp)
              
              Position(1) = bobvx(1,1) + rx * xx
              Position(2) = bobvx(1,2) + ry * yy
              Position(3) = bobvx(1,3) + rz * zz
              !
              !Find element containing point, this is the core of the voxelization
              !
              call elsest_host_element(&
                   ielse,relse,1_ip,meshe(ndivi),Position,ielem,&
                   shapp,deriv,coloc,dista)

              !call elsest(&
              !     2_ip,1_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),&
              !     lnods,ltype,ltopo,coord,Position,relse,&
              !     ielem,shapp,deriv,coloc,dummi)

              if( ielem > 0 ) then ! Something found
                 nvoxl=nvoxl+1
                 if (nvoxl>OldSize) then !There are more voxels than we guessed, shouldn't happen
                    ! We must reallocate memory to cope with more than expected voxels
                    allocate(idxvx(OldSize),elevx(OldSize))
                    idxvx = TempIndexList
                    elevx = TempElemList
                    deallocate(TempElemList,TempIndexList)
                    allocate(TempIndexList(2*OldSize),TempElemList(2*OldSize))
                    TempIndexList(1:Oldsize)=idxvx
                    TempElemList(1:Oldsize)=elevx
                    deallocate(idxvx,elevx)
                    OldSize = size(TempElemList)
                 endif

                 TempElemList(nvoxl)  = ielem !Which element contains the voxel
                 TempIndexList(nvoxl) = ix + (iy-1)*resvx(1) + (iz-1)*resvx(1)*resvx(2) !Index for the voxel

              endif
           enddo
        enddo
     enddo
  end if

  if(nvoxl>0) then !We transfer list of elements and index list to their final place
     allocate (elevx(nvoxl),idxvx(nvoxl))
     elevx = TempElemList(1:nvoxl) 
     idxvx = TempIndexList(1:nvoxl)
  else
     allocate (elevx(1),idxvx(1))   !! CHECK IF THIS IS CORRECT IN SEQUENTIAL MODE
     elevx=0
     idxvx=0
  end if
  
  if( INOTMASTER ) deallocate (TempElemList,TempIndexList) !Get rid of temps

  !----------------------------------------------------------------------
  !
  ! Gather: Master gathers the size to be sent by each slave
  !
  !----------------------------------------------------------------------

  if (ISLAVE) then
     !Transfer number of voxels found inside domain
     dummp(1) =  nvoxl  !We put info in an array
     paris    => dummp  !Send array
     parig    => nul1i  !Receive array is null for slaves
     npasi    =  1      !How many integers to send
  else if (IMASTER) then
     !Receive voxel number from slaves
     call memgen(1_ip,npart+1,0_ip)
     dummp(1) =  0  !Master sends nvox=0
     paris    => dummp !Send array
     parig    => gisca !Received data will be stored here
     npasi    =  1     
     npari    =  1
  end if

  call par_lagran(1_ip)       !This is a gather operation, each slave sends one int

  !
  ! Now the master knows how many voxels each slave has, it needs to receive
  ! info about WHICH voxels. This will be done with GatherV
  !
  if( IMASTER ) then
     nslvx = gisca
     !
     ! Construct displacement list
     !
     dplvx(1) = 0
     dplvx(2) = 0
     do is=3,npart+1
        dplvx(is) = dplvx(is-1)+nslvx(is-1)  ! equiv to Cumulative number of voxels per slave
     enddo
     call memgen(3_ip,npart+1,0_ip) !Deallocate memory
     !
     !Total number of voxels detected by slaves can be greater than actual number of voxels
     !
     totvx=dplvx(npart+1)+nslvx(npart+1)  
  else ! We are sequential
     ! I am not sure if the following will work in all cases
     totvx = product(resvx)
  end if

  !----------------------------------------------------------------------
  !
  ! Gatherv: Master gathers the voxels index list and put them into PARIN
  !
  !----------------------------------------------------------------------

  if (ISLAVE) then
     paris    => idxvx ! Send buffer contains the index of the voxels that each slave has
     npasi    =  nvoxl ! nvoxl can be zero from domains that did not detect voxels
     parin    => nul1i ! Not important for slaves
     parig    => nul1i
     pari1    => nul1i
  else if (IMASTER) then
     call memgen(1_ip,totvx,0_ip)
     allocate(gisca(totvx))  ! for some reason gisca is not allocated by memgen
     paris    => idxvx ! Master still needs to send something
     npasi    =  0     ! We know that master has nvoxl = 0
     parin    => gisca ! Receive buffer
     parig    => nslvx ! receive count array
     pari1    => dplvx ! Displacement list
  end if

  call par_lagran(4_ip) ! GatherV with integers

  !----------------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !----------------------------------------------------------------------

  if (IMASTER) then
     deallocate(idxvx)
     allocate(idxvx(totvx))
     idxvx=gisca 
     call memgen(3_ip,totvx,0_ip)
  endif

end subroutine prevox
