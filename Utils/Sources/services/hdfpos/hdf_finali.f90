subroutine hdf_finali()
  use def_kintyp
  use def_hdfpos
  use def_master  
  use mod_parall, only : PAR_INTEGER,PAR_COMM_MY_CODE
  implicit none

#ifdef MPI_OFF
#else
  include 'mpif.h'
#endif

  integer(4) :: ierr

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

#ifdef MPI_OFF
  call h5close_f( ierr )
#else
  if (hdf5_comm /= MPI_COMM_NULL) then

    !
    ! Close property lists.
    !
    call h5pclose_f( hdf5_faplist, ierr )
    call h5pclose_f( hdf5_dxplist, ierr )
    call h5pclose_f( hdf5_dcplist, ierr )

    !
    ! Close HDF5 library.
    !
    call h5close_f( ierr )

    ! Deallocate group and communicator for HDF5
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'deallocating comm and group:', hdf5_comm, hdf5_group
#endif
    call MPI_COMM_FREE( hdf5_comm, ierr )
    call MPI_GROUP_FREE( hdf5_group, ierr )

    ! Free MPI_Info
    call MPI_INFO_FREE( hdf5_info, ierr )
  endif
#endif

  !
  ! Deallocate HDF5 arrays
  !
  deallocate(hdf5_npoins)
  deallocate(hdf5_nelems)

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine hdf_finali

