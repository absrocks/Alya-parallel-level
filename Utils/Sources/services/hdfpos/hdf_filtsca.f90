subroutine hdf_filtsca
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


  integer(4)                :: i, j, ierr, rank = 1,tot
  integer(hid_t)            :: filespace
  integer(hid_t)            :: dset_id     ! Dataset identifier
  integer(hid_t)            :: memspace    ! Memspace identifier
  integer(hsize_t), pointer :: offset(:)
  integer(hsize_t)          :: dims(1)
  integer(ip),      pointer :: sizee(:)

  real(rp),         pointer :: gesca_filter(:)
  integer(ip)               :: ipoin
  character(7)              :: wopofilt

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  !
  ! Create the data space for the dataset.
  !
  dims(1) = npoin_total_filt
  write(wopofilt,'(a,i1)') 'FILTER',kfl_filte

  call h5screate_simple_f( rank, dims, filespace, ierr )

  !
  !  Create open a group 
  !
  if (.not.hdf5_group_filter) then
     call h5gcreate_f(hdf5_fileid,wopofilt, group_id, ierr)
     hdf5_group_filter = .true.
  end if

#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_filter: h5screate_simple_f filespace dims:', dims
#endif
  !
  ! Create the dataset with default properties.
  !
#ifdef R4
  call h5dcreate_f( group_id, wopos_hdf(1), H5T_NATIVE_REAL, filespace, &
                    dset_id, ierr, dcpl_id = hdf5_dcplist )
#else
  call h5dcreate_f( group_id, wopos_hdf(1), H5T_NATIVE_DOUBLE, filespace, &
                    dset_id, ierr, dcpl_id = hdf5_dcplist )
#endif
  call h5sclose_f( filespace, ierr )
  !
  ! Each process defines dataset in memory and writes it to the hyperslab
  ! in the file. 
  !
  dims(1) = npoin_filt
  !
  ! new array to put the filtered data (gesca_hdf>gesca_filter)
  !
  allocate(gesca_filter(npoin_filt))
  i=1_ip
  do ipoin=1,npoin
     if(gefil(ipoin)/=0) then
        gesca_filter(i)=gesca_hdf(ipoin)
        i=i+1
     end if
  end do

  call h5screate_simple_f( rank, dims, memspace, ierr )

#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_filter: h5screate_simple_f memspace dims:', dims
#endif

  ! 
  ! Select hyperslab in the file.
  !
  allocate(offset(hdf5_size))
  do i= 1, hdf5_size
     offset(i) = hdf5_filter(i)
  enddo

  call h5dget_space_f( dset_id, filespace, ierr )
  call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, offset(hdf5_rank+1:), dims, ierr )


#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_filtsca: offset(hdf5_rank+1)', offset(hdf5_rank+1)
  write(*,*) "hdf5_size",hdf5_size
  write(*,*) "hdf5_filter",hdf5_filter
  write(*,*) "kfl_paral, npoin_filt",kfl_paral,npoin_filt
#endif

  deallocate(offset)
  deallocate(hdf5_filter)

  !
  ! Write the dataset collectively. 
  !
#ifdef R4
  call h5dwrite_f( dset_id, H5T_NATIVE_REAL, gesca_filter, dims, ierr, &
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#else
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, gesca_filter, dims, ierr, &
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#endif
#ifdef HDF5_DEBUG
  do i=1, npoin_filt
     !write(*,*) kfl_paral, 'hdf_filter gevec_hdf:', gevec_hdf(1:SIZE(gevec_hdf, DIM=1), i)
     write(*,*) kfl_paral, 'hdf_filter gesca_filter:', gesca_filter(i)
  enddo
#endif

  !
  ! Close the group
  !
  !call h5gclose_f(group_id, ierr)
  !
  ! Close dataspaces.
  !
  call h5sclose_f( filespace, ierr )
  call h5sclose_f( memspace, ierr )

  !
  ! Close the dataset.
  !
  call h5dclose_f( dset_id, ierr )

  !
  !
  !
  deallocate(gesca_filter)


#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine hdf_filtsca
