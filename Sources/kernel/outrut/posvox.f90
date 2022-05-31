subroutine posvox(wopos)
  !-----------------------------------------------------------------------
  !****
  ! NAME
  !    posvox
  ! DESCRIPTION
  !    This routine saves data in a predefined box into voxel format.
  !    Each slave loops through the voxels contained in its domain, computing
  !    the value of the variable at the voxel point.
  !    Afterwards, slaves send this data (and only the data) to the master, who
  !    puts it in its correct order using the info from preprocessing (prevox.f90)
  !    and outputs it to file.
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_postpr
  use mod_elmgeo, only  :  elmgeo_natural_coordinates
  use mod_opfpos, only  :  opfposvx
  implicit none
  character(5), intent(in)  :: wopos(3)
  integer(ip)               :: ivoxl,ipoin,ielem,inode,pnode,pelty
  integer(ip)               :: idime,ix,iy,iz,ifoun,index
  real(4),      allocatable :: VoxelValues(:,:,:)
  real(rp)                  :: Position(3),lmini,lmaxi,gmini,gmaxi
  real(rp)                  :: shapp(mnode),deriv(mnode*3),coloc(3)
  real(rp)                  :: elcod(ndime,mnode),value,valu2(3)
  real(rp),     pointer     :: gtemp(:)
  integer(ip)               :: triple(3),idx
  integer(4)                :: resv4(3)
  !
  ! Don't know how to voxelize in 1 or 2 dims yet
  !
  if( ndime < 3 ) then
     return
  endif
  !
  ! Slaves or master without voxels still need to send an array of size one
  !
  if( nvoxl > 0 ) then
     allocate( vxdat(nvoxl) )
  else
     allocate( vxdat(1) )
     vxdat(1) = 0.0_rp
  endif
  !
  ! Automatic detection of maximum and minimum
  !
  gmini =  1.0e12_rp
  gmaxi = -1.0e12_rp
  lmini = -relse(1)
  lmaxi =  1.0_rp + relse(1)

  !----------------------------------------------------------------------
  !
  ! Interpolate values at my voxels
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then

     do ivoxl = 1,nvoxl
        !
        ! Loop over my voxels
        !
        Position = voxpos ( idxvx(ivoxl) )   ! Position of voxel
        ielem = elevx (ivoxl)                ! Element containing the voxel
        pelty = abs(ltype(ielem))            ! Prepare for calling elsest
        pnode = nnode(pelty)
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           do idime = 1,ndime
              elcod(idime,inode) = coord(idime,ipoin)
           end do
        end do
        !
        ! Interrogate elsest about the shape function at the given position
        !
        call elmgeo_natural_coordinates(    &
             ndime,pelty,pnode,elcod,shapp, &
             deriv,Position,coloc,ifoun)
        if( ifoun <= 0 ) then
           call runend('SOMEHOW A VOXEL PREVIOUSLY FOUND IN AN ELEMENT IS NOT THERE ANYMORE...')
        end if
        !
        ! Select which variable to save and compute it using the shape function
        !
        value = 0.0_rp
        if( wopos(2) == 'SCALA' ) then
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              value = value + shapp(inode) * gesca(ipoin)
           end do
        else if( wopos(2) == 'VECTO' ) then
           valu2(1) = 0.0_rp
           valu2(2) = 0.0_rp
           valu2(3) = 0.0_rp
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 valu2(idime) = valu2(idime) + shapp(inode) * gevec(idime,ipoin)
              end do
           end do
           do idime = 1,ndime
              value = value + valu2(idime) * valu2(idime)
           end do
           value = sqrt(value)
        else
           call runend('FOR NOW VOXELS CAN ONLY SAVE SCALARS OR VECTOR MODULES.')
        end if
        gmini = min(gmini,value)
        gmaxi = max(gmaxi,value)
        vxdat(ivoxl) = value
     end do

  end if
  !
  ! Compute global Min and max values
  !
  call pararr('MIN',0_ip,1_ip,gmini)
  call pararr('MAX',0_ip,1_ip,gmaxi)

  !----------------------------------------------------------------------
  !
  ! Gatherv voxel results
  !
  !----------------------------------------------------------------------

  if( ISLAVE ) then
     parrs => vxdat  ! Send buffer
     npasr =  nvoxl  ! how much data to send
     parre => nul1r
     parig => nul1i
     pari1 => nul1i
  else if( IMASTER ) then
     allocate( gtemp(totvx) )
     parrs => vxdat  ! Needed in any case
     npasr =  0      ! We know that master has nvoxl = 0
     parre => gtemp  ! Receive data will go here
     parig => nslvx  ! Receive count array
     pari1 => dplvx  ! Displacement list
  end if

  call par_lagran(2_ip)  ! GatherV with reals

  if( IMASTER ) then
     deallocate(vxdat)
     allocate(vxdat(totvx))
     do ix = 1,totvx
        vxdat(ix) = gtemp(ix)
     end do
     deallocate(gtemp)
  endif

  if( INOTSLAVE ) then

     !----------------------------------------------------------------------
     !
     ! Reconstruct voxel data tensor
     !
     !----------------------------------------------------------------------

     allocate(VoxelValues(resvx(1),resvx(2),resvx(3)))
     do iz= 1, resvx(3)
        do iy = 1, resvx(2)
           do  ix = 1, resvx(1)
              VoxelValues(ix,iy,iz) = 0.0_4
           end do
        end do
     end do
     do ivoxl = 1,totvx
        idx    = idxvx(ivoxl)  ! Index of the voxel
        triple = idx2tr(idx)   ! x,y,z values for the index
        ! We shift the value so that it is always positive -- this is important so that each
        ! individual voxel file can be read directly, but for animations one must pass
        ! by alya2blender.py to equalize scale among frames
        VoxelValues (triple(1),triple(2),triple(3)) = real(vxdat(ivoxl)-gmini,4_4)
     end do

     !----------------------------------------------------------------------
     !
     ! Output to binary file, Blender format
     !
     !----------------------------------------------------------------------

     wopos_pos(1) = wopos(1)
     call opfposvx()
     !
     ! Write header just for post-checking, it will be stripped by postprocessing if
     ! many frames will be stiched together.
     ! Precision must be 4 bytes for both integer and reals
     !
     resv4(1) = int(resvx(1),4)
     resv4(2) = int(resvx(2),4)
     resv4(3) = int(resvx(3),4)

     write (lun_posvx) resv4(1)
     write (lun_posvx) resv4(2)
     write (lun_posvx) resv4(3)
     write (lun_posvx) 1_4    ! This indicates just one frame, alya2blender.py uses it to check endiannes

     ! Save data
     do iz= 1, resvx(3)
        do iy = 1, resvx(2)
           do  ix = 1, resvx(1)
              write (lun_posvx) VoxelValues(ix,iy,iz)
           end do
        end do
     end do

     ! We do a little trick here, we save the true min and max values at the end of the file
     ! where Blender does not look, just in case we need to reconstruct the data
     ! or create a uniform scale among different frames.
     write (lun_posvx) real(gmini,4)
     write (lun_posvx) real(gmaxi,4)

     close(unit=lun_posvx)
     deallocate( VoxelValues )

  end if
  !
  ! Deallocate memory
  !
  deallocate( vxdat )

end subroutine posvox

