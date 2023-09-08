subroutine hdf_openfi(itask)
  !-----------------------------------------------------------------------
  !****f* hdfpos/hdf_openfi
  ! NAME
  !    hdf_openfi
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Sends and saves postprocess file name
  !    ITASK = 2 ... Composes the file name when using HDF format
  ! USED BY
  !    Hdfpos
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_hdfpos
  use def_domain
  use def_postpr
  use mod_parall, only : PAR_INTEGER,PAR_COMM_MY_CODE
  !
  !  use h5fddsm
  !
  implicit none
#ifdef MPI_OFF
#else
  include 'mpif.h'
#endif
  integer(ip),    intent(in) :: itask
  character(150), save       :: filsa

  integer(ip)                :: kfl_ptask_old
  integer(4)                 :: ierr

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  select case ( itask )

  case ( 1_ip )

#ifndef MPI_OFF
     !-------------------------------------------------------------------
     !
     ! Send and save postprocess file name
     !
     !-------------------------------------------------------------------
     !
     ! Master sends file name to slaves
     !
     kfl_ptask_old = kfl_ptask
     kfl_ptask = 1
     strch = 'hdf_openfi'
     do parii = 1,2 
        npari = 0
        nparr = 0
        nparc = 0
        nparc = nparc + 150
        if( parii == 2 .and. IMASTER ) parch(1:150)     = fil_postp(1:150)
        if( parii == 2 .and. ISLAVE  ) fil_postp(1:150) = parch(1:150) 
        if( parii == 1 ) then
           if( ISLAVE ) then
              call par_broadc()
           end if
        end if
     end do
     if( IMASTER ) call par_broadc()
#endif
     !
     ! Save file name
     !
     filsa = trim(fil_postp)
     kfl_ptask = kfl_ptask_old

  case ( 2_ip )

     !-------------------------------------------------------------------
     !
     ! Compose postprocess file name for output data
     !
     !-------------------------------------------------------------------
     !
     ! Time step identifier 
     !
     write(nunam_pos,'(i7)') ittim         
     if( ittim < 10 ) then
        write(nunam_pos,'(a,i1)') '000000',ittim
     else if( ittim < 100 ) then
        write(nunam_pos,'(a,i2)') '00000',ittim
     else if( ittim < 1000 ) then
        write(nunam_pos,'(a,i3)') '0000',ittim
     else if( ittim < 10000 ) then
        write(nunam_pos,'(a,i4)') '000',ittim
     else if( ittim < 100000 ) then
        write(nunam_pos,'(a,i5)') '00',ittim
     else if( ittim < 1000000 ) then
        write(nunam_pos,'(a,i6)') '0',ittim
     end if
     !
     ! Compose file name: example-1234567.h5
     !
     fil_postp = trim(filsa)//'-'//adjustl(trim(nunam_pos))//'.h5'


     !
     ! Open the file collectively
     !

     !
     ! Create the file collectively.
     ! 
     call h5fcreate_f( fil_postp, H5F_ACC_TRUNC_F, hdf5_fileid, ierr, access_prp=hdf5_faplist )
#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_openfi: opening h5 file:', hdf5_fileid
#endif

  case ( 3_ip )

     !-------------------------------------------------------------------
     !
     ! Compose postprocess file name for mesh
     !
     !-------------------------------------------------------------------
     ! Compose file name: example-mesh.h5
     !
     fil_postp = trim(filsa)//'-mesh.h5'

  case ( 4_ip )

     !-------------------------------------------------------------------
     !
     ! Compose restart file name for output data
     !
     !-------------------------------------------------------------------
     !
     ! Time step identifier 
     !
     write(nunam_pos,'(i7)') ittim         
     if( ittim < 10 ) then
        write(nunam_pos,'(a,i1)') '000000',ittim
     else if( ittim < 100 ) then
        write(nunam_pos,'(a,i2)') '00000',ittim
     else if( ittim < 1000 ) then
        write(nunam_pos,'(a,i3)') '0000',ittim
     else if( ittim < 10000 ) then
        write(nunam_pos,'(a,i4)') '000',ittim
     else if( ittim < 100000 ) then
        write(nunam_pos,'(a,i5)') '00',ittim
     else if( ittim < 1000000 ) then
        write(nunam_pos,'(a,i6)') '0',ittim
     end if
     !
     ! Compose file name: example-1234567.h5
     !
     if (kfl_rsfil == 1) then
        fil_rstar = trim(filsa)//'-'//adjustl(trim(nunam_pos))//'.rst.h5'
     else 
        fil_rstar= trim(filsa)//'.rst.h5'
     end if

     !
     ! Open the file collectively
     !

     !
     ! Create the file collectively.
     ! 
     call h5fcreate_f( fil_rstar, H5F_ACC_TRUNC_F, hdf5_fileid, ierr, access_prp=hdf5_faplist )
#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_openfi: opening h5 restar file:', hdf5_fileid
#endif

  case ( 5_ip )

     !-------------------------------------------------------------------
     !
     ! Recompose restart file name 
     !
     !-------------------------------------------------------------------
     ! 
     !
     fil_rstar= adjustl(trim(filsa))//'.rst.h5'


  end select

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine hdf_openfi

