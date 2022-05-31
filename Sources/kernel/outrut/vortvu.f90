!------------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @file    vortvu.f90
!> @date    31/08/2012
!> @author  Hadrien Le Grand (Con) Calmet
!> @brief   pipi pipi
!> @details  pipi pipi
!> @}
!------------------------------------------------------------------------
subroutine vortvu

  use def_kintyp
  use def_master
  use def_domain
  use def_kermod

  implicit none


  integer(ip)             :: ii,is,totnp
  !real(rp)                :: 

  character(8)            :: chtim 
  integer(ip), target     :: dummp(1)
  integer(ip), pointer    :: icount(:),dplnp(:)
  real(rp),    pointer    :: coord_jpt(:),npdat(:)


  if (ittim == 0) return
  

  if (INOTSLAVE) then ! Master allocates memory for indexes of displacements and number of points
     allocate (icount(npart+1),dplnp(npart+1)) 
     !
     ! file for vu , writing the header of .vu
     !
     if(ittim<10) then
        write(chtim,'(a,i1)') '0000000',ittim
     else if(ittim<100) then
        write(chtim,'(a,i2)') '000000',ittim
     else if(ittim<1000) then
        write(chtim,'(a,i3)') '00000',ittim
     else if(ittim<10000) then
        write(chtim,'(a,i4)') '0000',ittim
     else if(ittim<100000) then
        write(chtim,'(a,i5)') '000',ittim
     else if(ittim<1000000) then
        write(chtim,'(a,i6)') '00',ittim
     else if(ittim<10000000) then
        write(chtim,'(a,i7)') '0',ittim
     end if

     open(unit=51,file=trim(chtim)//'.res.vu',status='unknown')
     write(51,*)"CHAMP Coo( ) ="
     write(51,*)"{"
     write(51,*)"// DonnÃ©es x, y, z" 

  end if

  if(ISLAVE) then
     !
     !allocation the points array with the right size for each proc 
     !and npdat equal but for the mpi communication
     !
     allocate (npdat(nvort*ndime))

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
  end if
!!!!!!!mpi_gather!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call par_lagran(1_ip)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  else ! We are sequential
     totnp = nvort
  end if
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
     allocate( coord_jpt(totnp))
     parrs => nul1r
     npasr =  0      ! We know that master has nvoxl = 0
     parre => coord_jpt
     parig => icount  ! Receive count array
     pari1 => dplnp  ! Displacement list
  end if
!!!!!!!!!MPI_GATHERV!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call par_lagran(2_ip)    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (IMASTER) then
     do ii=1,totnp,3
        write(51,*)coord_jpt(ii),coord_jpt(ii+1),coord_jpt(ii+2)
     end do

     write(51,*)"};"
     write(51,*)
     write(51,*)"MAILLAGE MonMaillage( ) ="
     write(51,*)"{"
     write(51,*)"   ZONE Zone1( LagrPoint01, Coo%3,",totnp/3,");"
     write(51,*)"};"      
     close(51) 

     deallocate(coord_jpt)
     deallocate (icount) 
     deallocate (dplnp)
  end if


end subroutine vortvu
