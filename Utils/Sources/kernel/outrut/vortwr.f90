!------------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @file    vortwr.f90
!> @date    07/10/2014
!> @author  Hadrien Calmet
!> @brief   Alya output core vortex
!> @details Write in .csv file the core of the vortex 
!> @}
!------------------------------------------------------------------------
subroutine vortwr

  use def_kintyp
  use def_master
  use def_domain
  use def_kermod
  use mod_memory
  use mod_communications, only : PAR_GATHER
  use mod_communications, only : PAR_GATHERV
  implicit none 

  integer(ip)             :: ii,is,totnp,istat
  character(8)            :: chtim 
  integer(ip), target     :: dummp(1)
  integer(ip), pointer    :: icount(:),dplnp(:)
  real(rp),    pointer    :: coord_jpt(:),npdat(:)
  integer(ip)             :: filevor = 888

  character(150)          :: filso
  character(150)          :: nunam_pos1

  !
  ! no data writing at T=0
  !
  if (ittim == 0) return
  !
  ! Master allocates memory for indexes of displacements and number of points
  !
  if (INOTSLAVE) then 
     !
     nullify (icount)
     nullify (dplnp) 
     !
     call memory_alloca(memor_dom,'ICOUNT','vortwr',icount,npart+1_ip)
     call memory_alloca(memor_dom,'DPLNP', 'vortwr',dplnp, npart+1_ip)
     !
     !setting of the output name
     !
     write(nunam_pos1,'(i7)') ittim         
     if( ittim < 10 ) then
        write(nunam_pos1,'(a,i1)') '000000',ittim
     else if( ittim < 100 ) then
        write(nunam_pos1,'(a,i2)') '00000',ittim
     else if( ittim < 1000 ) then
        write(nunam_pos1,'(a,i3)') '0000',ittim
     else if( ittim < 10000 ) then
        write(nunam_pos1,'(a,i4)') '000',ittim
     else if( ittim < 100000 ) then
        write(nunam_pos1,'(a,i5)') '00',ittim
     else if( ittim < 1000000 ) then
        write(nunam_pos1,'(a,i6)') '0',ittim
     end if
     filso = trim(title)//'-VORTX-'//trim(nunam_pos1)
     !
     ! open the output name (VTK or NOT VTK)
     !
     if ( kfl_outfo == 40 .OR. kfl_outfo == 41) then
        open (unit=filevor,file=trim(title)//'_vtk/'//trim(filso)//'.csv',status='unknown')
        write(filevor,'(a)') ' xcoord, ycoord, zcoord '
     else
        open (unit=filevor,file=trim(filso)//'.csv',status='unknown')
        write(filevor,'(a)') ' xcoord, ycoord, zcoord '
     endif
  end if

  if(ISLAVE) then
     !
     !allocation the points array with the right size for each proc 
     !and npdat equal but for the mpi communication
     !
     nullify (npdat)
     call memory_alloca(memor_dom,'NPDAT','vortwr',npdat,nvort*ndime)
     !
     !
     do ii=1,nvort
        npdat(((ii-1)*3)+1)  =  gevec(1,ii)
        npdat(((ii-1)*3)+2)  =  gevec(2,ii)
        npdat(((ii-1)*3)+3)  =  gevec(3,ii)
     end do
     !
     !Transfer number of point found for each proc
     !      
     dummp(1) =  nvort
     paris    => dummp 
     parig    => nul1i
     npasi    =  1     
     !
     !master receive number of point from slaves
     !
  else if (IMASTER) then
     call memgen(1_ip,npart+1,0_ip)
     dummp(1) =  0  !Master sends npts=0
     paris    => dummp
     parig    => gisca
     npasi    =  1
     npari    =  1
  else if (ISEQUEN) then
     !
     ! 
     !
     do ii=1,nvort
        write(filevor,13) gevec(1,ii),gevec(2,ii),gevec(3,ii)
     end do
     close(filevor)
  end if
  !
  ! mpi_gather
  !
  call Parall(705_ip)


!  call PAR_GATHER(nvort,gisca,'IN MY CODE') 
!  call PAR_GATHERV(lneig,lneig_gat,nneig4_gat,'IN MY CODE') ! REGARDER kernel/parall/par_output_partition.f90


  !
  !Construct displacement list
  !
  if( IMASTER ) then
     icount = gisca*ndime
     dplnp(1) = 0
     dplnp(2) = 0
     do is=3,npart+1
        dplnp(is) = dplnp(is-1)+icount(is-1)
     enddo
     call memgen(3_ip,npart+1,0_ip)
     !
     !Total number of points detected by slaves 
     !
     totnp=dplnp(npart+1)+icount(npart+1)
     ! 
  else ! We are sequential
     totnp = nvort
  end if
  !
  !
  !
  
  !
  !Receive points index list from slaves
  !
  if( ISLAVE ) then
     parrs => npdat
     npasr =  nvort*ndime  
     parre => nul1r
     parig => nul1i
     pari1 => nul1i
  else if( IMASTER ) then
     !
     nullify (coord_jpt)
     call memory_alloca(memor_dom,'COORD_JPT','vortwr',coord_jpt,totnp)
     !
     parrs => nul1r
     npasr =  0         ! We know that master has nvoxl = 0
     parre => coord_jpt
     parig => icount    ! Receive count array
     pari1 => dplnp     ! Displacement list
  end if
  !
  ! mpi_gatherv
  !
  call Parall(706_ip)
  ! 
  ! master writes the resultrs
  !
  if( IMASTER ) then
     do ii = 1,totnp,3
        write(filevor,13) coord_jpt(ii),coord_jpt(ii+1),coord_jpt(ii+2)
     end do
     close(filevor)
     call memory_deallo(memor_dom,'COORD_JPT','vortwr',coord_jpt)
  endif
  !
  ! deallocation for master or sequencial
  !
  if(INOTSLAVE) then
     call memory_deallo(memor_dom,'ICOUNT','vortwr',icount)
     call memory_deallo(memor_dom,'DPLNP' ,'vortwr',dplnp)
  end if


13 format (e12.6,',',e12.6,',',e12.6)

end subroutine vortwr
