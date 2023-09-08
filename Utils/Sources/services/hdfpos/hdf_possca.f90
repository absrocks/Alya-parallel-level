subroutine hdf_possca()
  use def_kintyp
  use def_master
  use def_hdfpos
  use def_parall
  use def_domain
  use mod_parall, only : PAR_INTEGER,PAR_COMM_MY_CODE
  implicit none
#ifdef MPI_OFF
#else
  include 'mpif.h'
#endif

  integer(4)                :: i, ierr, rank = 1
  integer(hid_t)            :: filespace   ! Filespace identifier
  integer(hid_t)            :: dset_id     ! Dataset identifier
  integer(hid_t)            :: memspace    ! Memspace identifier
  integer(hsize_t)          :: dims(1)



#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  !
  ! Create the data space for the dataset. 
  !
  dims(1) = npoin_total
  call h5screate_simple_f( rank, dims, filespace, ierr )
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_possca: h5screate_simple_f filespace dims: ', dims, 'ierr: ', ierr
#endif
  if (ierr .lt. 0) then
    return
  endif

  !
  ! Create the dataset with default properties.
  !
#ifdef R4
  call h5dcreate_f( hdf5_fileid, wopos_hdf(1), H5T_NATIVE_REAL, filespace, &
                    dset_id, ierr, dcpl_id = hdf5_dcplist )
#else
  call h5dcreate_f( hdf5_fileid, wopos_hdf(1), H5T_NATIVE_DOUBLE, filespace, &
                    dset_id, ierr, dcpl_id = hdf5_dcplist )
#endif
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_possca: h5dcreate_f wopos_hdf: ', wopos_hdf, 'ierr: ', ierr
#endif
  if (ierr .lt. 0) then
    return
  endif

  call h5sclose_f( filespace, ierr )

  !
  ! Each process defines dataset in memory and writes it to the hyperslab
  ! in the file. 
  !
  dims(1) = npoin
  call h5screate_simple_f( rank, dims, memspace, ierr )
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_possca: h5screate_simple_f memspace dims: ', dims, 'ierr: ', ierr
#endif
  if (ierr .lt. 0) then
    return
  endif

  ! 
  ! Select hyperslab in the file.
  !
  call h5dget_space_f( dset_id, filespace, ierr )
  if (ierr .lt. 0) then
    return
  endif
  call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, hdf5_npoins(hdf5_rank+1:), dims, ierr )
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_possca: h5sselect_hyperslab_f hdf5_npoins(hdf5_rank):', &
                        hdf5_npoins(hdf5_rank+1), 'ierr: ', ierr
#endif
  if (ierr .lt. 0) then
    return
  endif
 
  !
  ! Write the dataset collectively. 
  !
#ifdef R4
  call h5dwrite_f( dset_id, H5T_NATIVE_REAL, gesca_hdf, dims, ierr, &
                   file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#else
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, gesca_hdf, dims, ierr, &
                   file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#endif
  if (ierr .lt. 0) then
    return
  endif

  !
  ! Close dataspaces.
  !
  call h5sclose_f( filespace, ierr )
  call h5sclose_f( memspace, ierr )

  !
  ! Close the dataset.
  !
  call h5dclose_f( dset_id, ierr )

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine hdf_possca
