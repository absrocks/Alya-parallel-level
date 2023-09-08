SUBROUTINE Hdfpos(itask)
  use def_kintyp
  use def_master
  use def_domain
  use def_hdfpos
  use def_kermod 
  implicit none
#ifdef MPI_OFF
#else
  include 'mpif.h'
#endif

  integer(ip), intent(in) :: itask
  integer(ip)             :: i,j,ipoin
  real(rp)       :: time1,time2,time3
  integer,dimension(8) :: values
  real(rp)       :: dt1,dt2,dt3

  i=0_ip
  j=0_ip

#ifdef EVENT
  call mpitrace_user_function(1)
#endif
  
  select case ( itask )

  case ( 0_ip )
     !
     ! Read service
     !
     servi = ID_HDFPOS
     call reaser()

  case ( 1_ip )
    call date_and_time(VALUES=values)
    dt1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
    call cpu_time(time1)
     !
     ! Postprocess a scalar   
     !
     !   restart case
     !
     call hdf_restar()
     if ( kfl_reawr == 1 ) return
     if ( kfl_reawr == 2 ) return
     !
     !   get filter
     !
     call fildef(3_ip)
     !
     ! all gather npoin filtered for each slave
     !
     if (kfl_filte/=0.and.ISLAVE)then
        npoin_filt=0_ip
        do ipoin=1,npoin
           if(gefil(ipoin)/=0) npoin_filt=npoin_filt+1_ip
        end do
        call hdf_filter(0_ip)
     end if
     !
     !MPI_REDUCE : add npoin_filt for each slave and put in npoin_total_filt in master  
     !
     if (kfl_filte/=0) call hdf_filter(1_ip)
     !
     !
     !
     if (.not. hdf5_opened) then
        !
        ! Open H5 timestep output file
        !
        hdf5_opened = .true.
        if (ISLAVE) then
           call hdf_openfi(2_ip)
        endif
     endif
     !
     ! writing mesh xmf for each timestep
     !    
     if (IMASTER .and. i==0_ip .and. kfl_filte==0) then
        call hdf_wrtxmf(6_ip)
        i=1_ip
     else if (IMASTER .and. j==0_ip .and. kfl_filte/=0) then
        call hdf_wrtxmf(12_ip)
        j=1_ip 
     end if

     if (ISLAVE) then
        if (kfl_filte==0)then
           call hdf_possca()
        else
           call hdf_filtsca()
        end if
        if (.not. hdf5_mesh) then
           hdf5_mesh = .true.
           call hdf_posgeo(2_ip)
        end if
        !if (.not. hdf5_mesh_filter) then
        !   hdf5_mesh_filter = .true.
        !   write(*,*)'writing filter mesh h5'
        !end if

     else
        if (kfl_filte==0)then 
           call hdf_wrtxmf(2_ip)
        else
           call hdf_wrtxmf(13_ip)
        end if
     endif
     call date_and_time(VALUES=values)
     dt2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
     call cpu_time(time2)
     ! if (kfl_paral == 25) then
     !   write(*,*) 'TIME HDF5 1:' , time2-time1, ' ', dt2 - dt1
     ! end if
  case ( 2_ip )
     call date_and_time(VALUES=values)
     dt1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
     call cpu_time(time1)
     !
     ! Postprocess a vector
     !
     !   restart case
     call hdf_restar()
     if ( kfl_reawr == 1 ) return
     if ( kfl_reawr == 2 ) return
     !
     !   get filter
     call fildef(3_ip)
     !
     ! all gather npoin filtered for each slave
     !
     if (kfl_filte/=0.and.ISLAVE)then
        npoin_filt=0_ip
        do ipoin=1,npoin
           if(gefil(ipoin)/=0) npoin_filt=npoin_filt+1_ip
        end do
        call hdf_filter(0_ip)
     end if
     !
     !MPI_REDUCE : add npoin_filt for each slave and put in npoin_total_filt in master  
     !
     if (kfl_filte/=0) call hdf_filter(1_ip) 
     !
     !
     !
     if (.not. hdf5_opened) then
        !
        ! Open H5 timestep output file
        !
        hdf5_opened = .true.
        if (ISLAVE) then
           call hdf_openfi(2_ip)
        endif
     endif
     !
     ! writing mesh xmf for each timestep
     !    
     if (IMASTER .and. i==0_ip.and.kfl_filte==0) then
        call hdf_wrtxmf(6_ip)
        i=1_ip
     else if (IMASTER .and. j==0_ip.and.kfl_filte/=0) then
        call hdf_wrtxmf(12_ip)
        j=1_ip
     end if

     if (ISLAVE) then
        if (kfl_filte==0)then
           call hdf_posvec()
        else
           call hdf_filtvec()
        end if
     else
        if (kfl_filte==0)then 
           call hdf_wrtxmf(3_ip)
        else
           call hdf_wrtxmf(14_ip)
        end if
     endif
     call date_and_time(VALUES=values)
     dt2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
     call cpu_time(time2)
     ! if (kfl_paral == 25) then
     !   write(*,*) 'TIME HDF5 2:' , time2-time1, ' ', dt2-dt1
     ! end if
  case ( 3_ip )
     !call cpu_time(time1)
     !
     ! Postprocess a r3p 
     !
     !
     !   restart case
     !
     call hdf_restar()
     if ( kfl_reawr == 1 ) return
     if ( kfl_reawr == 2 ) return
     !
     ! no post process
     !

     !call cpu_time(time2)
     !if (kfl_paral == 25) then
     !  write(*,*) 'TIME HDF5 3:' , time2-time1
     !end if
  case ( 7_ip )
     call date_and_time(VALUES=values)
     dt1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
     call cpu_time(time1)
     !
     ! Close H5 timestep output file
     !
     if (hdf5_opened) then

        hdf5_opened       = .false.
        hdf5_mesh         = .false.
        hdf5_mesh_filter  = .false. 
        hdf5_group_filter = .false.
        hdf5_restar_write = .false.
        hdf5_restar_read  = .false.

        if (ISLAVE) then
           if (kfl_filte==0)then
              call hdf_closfi()
           else
              call hdf_filter(2_ip)
           end if
        else
           call hdf_wrtxmf(7_ip)
           !call hdf_wrtxmf(14_ip)????? backward 
        endif
     endif
     call date_and_time(VALUES=values)
     dt2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds
     call cpu_time(time2)
     ! if (kfl_paral == 25) then
     !   write(*,*) 'TIME HDF5 7:' , time2-time1, ' ', dt2 - dt1
     ! end if
  case ( 8_ip )
     !
     ! Allgather npoins information among HDF5 slaves
     ! 
     if( .not. PART_AND_WRITE() .and. kfl_outfo == 50 .and. IPARALL ) then
        if (ISLAVE) then
           call hdf_posgeo(0_ip)
        endif
     endif

     !
     ! Postprocess Coordinates and Connectivities
     !
     !          test.mesh.h5
     !                /
     !                |
     !       -------------------
     !       |                 |
     !     COORD             CONNE
     !
     call hdf_openfi(3_ip)
     if (ISLAVE) then
        call hdf_posgeo(1_ip)
     endif
     !
     ! XDMF writing the entete of mesh and filter
     !
     call hdf_wrtxmf(5_ip)
     !
     ! case of filter -> writing the entete
     !
     if ( kfl_filte /= 0 ) call hdf_wrtxmf(11_ip)
     !
     ! case of vortex core -> writing the entete
     ! 
     if ( kfl_vortx /= 0 ) call hdf_wrtxmf(8_ip)

  case ( 9_ip )
     !
     ! Master sends postprocess file name to slaves
     !
     call hdf_openfi(1_ip)
     !
     ! Initialize HDF
     ! - Create communicator
     ! - Initialize FORTRAN predefined data types
     !
     call hdf_initia()


  case ( 10_ip )
     !
     ! XDMF writing the end-file 
     !
     if (kfl_outfo == 50) then
        if (IMASTER) then
           call hdf_wrtxmf(4_ip)
           !
           ! XDMF writing the end-file of filer file
           !
           if ( kfl_filte /= 0 ) call hdf_wrtxmf(15_ip)
           !
           ! XDMF writing the end-file of vortex core
           !
           if ( kfl_vortx /= 0 ) call hdf_wrtxmf(10_ip)
           !
           ! XDMF writing the end-file of filer file
           !
           if (hdf5_mesh_filter) call hdf_wrtxmf(14_ip)
        end if
     end if
     !
     ! Finalize HDF
     !
     if (kfl_outfo == 50) then
        call hdf_finali()
     end if

  case ( 11_ip )
     !call cpu_time(time1)
     !
     !case of time = 0 no writing text neither hd5 data
     !
     if (ittim==0)then
        return
     endif
     !
     ! XDMF writing the entete of vortex core
     !
     !if (ittim == 1) then
     !   call hdf_wrtxmf(8_ip)
     !end if
     !
     ! Postprocess Coordinates of the vortex core
     ! Allgather nvort information among HDF5 slaves
     !
     if( .not. PART_AND_WRITE() .and. kfl_outfo == 50 .and. IPARALL ) then
        if (ISLAVE) then
           call hdf_posvor(0_ip)
           !write(*,*)'call hdf_posvor(0_ip)',kfl_paral
        endif
     endif
     !
     if (.not. hdf5_opened) then
        !
        ! Open H5 timestep output file
        !
        hdf5_opened = .true.
        if (ISLAVE) then
           call hdf_openfi(2_ip)
        endif
     endif
     !
     ! MPI_REDUCE (nvort_total -> master)
     !
     call hdf_posvor(2_ip)
     !write(*,*)'call hdf_posvor(2_ip)',kfl_paral
     !
     if (ISLAVE) then
        call hdf_posvor(1_ip)
        !write(*,*)'call hdf_posvor(1_ip)',kfl_paral
     else
        call hdf_wrtxmf(9_ip)
     endif
     !call cpu_time(time2)
     !if (kfl_paral == 25) then
     !  write(*,*) 'TIME HDF5 11:' , time2-time1
     !end if
  case( 12_ip )
     !
     ! Close restart file(call coming from /kernel/master/restart)
     !
     if (ISLAVE.and.hdf5_restar_write) call hdf_closfi()

  end select  

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine Hdfpos
