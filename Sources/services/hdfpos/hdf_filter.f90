subroutine hdf_filter(itask)
  use def_kintyp
  use def_master
  use def_hdfpos
  use def_parall
  use def_domain
  use def_kermod
  use mod_parall, only : PAR_INTEGER,PAR_COMM_MY_CODE
  implicit none
#ifdef MPI_OFF
#else
  include 'mpif.h'
#endif

  integer(ip), intent(in)   :: itask
  integer(4)                :: i, ierr
  integer(ip),      pointer :: sizee(:)
   
#ifdef EVENT
  call mpitrace_user_function(1)
#endif


  select case ( itask )

  case ( 0_ip )

     !
     ! Gather npoin filtered
     !
     allocate(sizee(hdf5_size))
     allocate(hdf5_filter(hdf5_size))
     npoin_total_filt=0_ip

#ifdef MPI_OFF
#else
     call MPI_ALLGATHER( npoin_filt, 1, PAR_INTEGER, &
          sizee, 1, PAR_INTEGER, hdf5_comm, ierr )
#endif

     hdf5_filter(1) = 0
     do i=2, hdf5_size
        hdf5_filter(i) = hdf5_filter(i-1) + sizee(i-1)
     enddo
     do i=1,hdf5_size
        npoin_total_filt = npoin_total_filt + sizee(i)
     end do

#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_filter: hdf5_filter', hdf5_filter
     write(*,*) kfl_paral, 'hdf_filter: npoin_total_filt', npoin_total_filt
#endif


     deallocate(sizee)

  case ( 1_ip )
     !
     ! MPI_REDUCE : add npoin_filt for each slave and put in npoin_total_filt in master
     !

#ifdef MPI_OFF
#else
     call MPI_REDUCE( npoin_filt, npoin_total_filt, 1, PAR_INTEGER, &
                      MPI_SUM, 0, PAR_COMM_MY_CODE, ierr )
#endif

#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_filter: npoin_total_filt', npoin_total_filt
#endif

  case ( 2_ip )
     !
     ! Close group H5 timestep file
     !   
     call h5gclose_f(group_id, ierr)
     
  end select


#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine hdf_filter
