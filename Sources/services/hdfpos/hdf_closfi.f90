subroutine hdf_closfi()
  !-----------------------------------------------------------------------
  !****f* hdfpos/hdf_closfi
  ! NAME
  !    hdf_openfi
  ! DESCRIPTION
  !    This routine:
  !    Close H5 file for the current timestep
  ! USED BY
  !    Hdfpos
  !***
  !-----------------------------------------------------------------------
  use def_hdfpos
  use def_parall
  use def_master
  implicit none
  integer(4)                :: ierr

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  !
  ! Close the file.
  !
  call h5fclose_f( hdf5_fileid, ierr )

#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_closfi: closing h5 file:', hdf5_fileid
#endif

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine hdf_closfi
