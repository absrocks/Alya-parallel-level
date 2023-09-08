subroutine hdf_posvec()
  use def_kintyp
  use def_master
  use def_hdfpos
  use def_parall
  use def_domain
  implicit none
#ifdef MPI_OFF
#else
  include 'mpif.h'
#endif

  integer(4)                :: i, j, ierr, rank = 2
  integer(hid_t)            :: filespace
  integer(hid_t)            :: dset_id     ! Dataset identifier
  integer(hid_t)            :: memspace    ! Memspace identifier
  integer(hsize_t), pointer :: offset(:,:)
  integer(hsize_t)          :: dims(2)
  real(rp)       :: time1pv,time2pv

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  !
  ! Create the data space for the dataset. 
  !
  dims(1) = SIZE(gevec_hdf, DIM=1)
  dims(2) = npoin_total
  call h5screate_simple_f( rank, dims, filespace, ierr )
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_posvec: h5screate_simple_f filespace dims:', dims
#endif
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
  call h5sclose_f( filespace, ierr )

  !
  ! Each process defines dataset in memory and writes it to the hyperslab
  ! in the file. 
  !
  dims(1) = SIZE(gevec_hdf, DIM=1)
  dims(2) = npoin
  call h5screate_simple_f( rank, dims, memspace, ierr )
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_posvec: h5screate_simple_f memspace dims:', dims
#endif

  !
  ! Select hyperslab in the file.
  !
  allocate(offset(2,hdf5_size))
  do i= 1, hdf5_size
    offset(1,i) = 0
    offset(2,i) = hdf5_npoins(i)
  enddo
  call h5dget_space_f( dset_id, filespace, ierr )
  call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, offset(1:2,hdf5_rank+1), dims, ierr )
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_posvec: offset(hdf5_rank+1)', offset(1:2,hdf5_rank+1)
#endif
  deallocate(offset)

  !
  ! Write the dataset collectively. 
  !
  !call cpu_time(time1pv)
#ifdef R4
  call h5dwrite_f( dset_id, H5T_NATIVE_REAL, gevec_hdf, dims, ierr, &
                   file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#else
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, gevec_hdf, dims, ierr, &
                   file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#endif
  !call cpu_time(time2pv)
  !if (kfl_paral == 25) then
  !   write(*,*) 'TIME HDF5 pv:' , time2pv-time1pv
  !end if
#ifdef HDF5_DEBUG
  do i=1, npoin
    write(*,*) kfl_paral, 'hdf_posvec: gevec_hdf:', gevec_hdf(1:SIZE(gevec_hdf, DIM=1), i)
  enddo
#endif

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

end subroutine hdf_posvec
