subroutine hdf_posr3p()
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

  integer(4)                :: i, ielem, ibuf, ierr, rank = 1
  integer(hid_t)            :: filespace
  integer(hid_t)            :: dset_id     ! Dataset identifier
  integer(hid_t)            :: memspace    ! Memspace identifier
  integer(hsize_t)          :: dims(1)
  integer(hsize_t), pointer :: offset(:)
  integer(hsize_t), pointer :: sizes(:)
  real(rp),         pointer :: buffer(:)

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  !
  ! Create the data space for the dataset. 
  !
  ! Calculate the required size for the buffer

  allocate(sizes(hdf5_size))
  dims(1) = 0
  do ielem = 1, nelem
     dims(1) = dims(1) + SIZE(ger3p_hdf(ielem)%a, DIM=1)*SIZE(ger3p_hdf(ielem)%a, DIM=2)
  enddo
  !write(*,*)'dims(1) =',dims(1)
  !    
  !    !!!!!!!!! case of post processing the vesgs with actual time and previus time !!!!!!!!
  !     
  !    do ielem = 1, nelem 
  !    dims(1) = dims(1) + SIZE(ger3p_hdf(ielem)%a)
  !    enddo   
  !
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, dims(1), SIZE(ger3p_hdf(ielem)%a), SIZE(ger3p_hdf(ielem)%a, DIM=1), &
       SIZE(ger3p_hdf(ielem)%a, DIM=2), SIZE(ger3p_hdf(ielem)%a, DIM=3)
#endif
  !
  ! Copy r3p data in a sequential buffer
  !
  allocate(buffer(dims(1))) ! In MN we have to add 1 
  ibuf = 1
  do ielem = 1, nelem    
     buffer(ibuf:ibuf + SIZE(ger3p_hdf(ielem)%a, DIM=1)*SIZE(ger3p_hdf(ielem)%a, DIM=2) - 1) & 
          = RESHAPE(ger3p_hdf(ielem)%a, (/SIZE(ger3p_hdf(ielem)%a, DIM=1)*SIZE(ger3p_hdf(ielem)%a, DIM=2)/))
     ibuf = ibuf + SIZE(ger3p_hdf(ielem)%a, DIM=1)*SIZE(ger3p_hdf(ielem)%a, DIM=2)          
  enddo
  !    
  !    !!!!!!!!! case of post processing the vesgs with actual time and previus time !!!!!!!!
  !
  !    buffer(ibuf:ibuf+SIZE(ger3p_hdf(ielem)%a)) = RESHAPE(ger3p_hdf(ielem)%a, &
  !                                                         (/SIZE(ger3p_hdf(ielem)%a)/))
  !    ibuf = ibuf + SIZE(ger3p_hdf(ielem)%a)
  !
  !
  ! Sanity check for memory copy
  if ((ibuf-1) .ne. dims(1)) then
     call runend('Error in hdf_posr3p copying r3p structure to sequential buffer')
  endif

#ifdef MPI_OFF
#else
  call MPI_ALLGATHER( dims, 1, MPI_INTEGER8, &
       sizes, 1, MPI_INTEGER8, hdf5_comm, ierr )
#endif

  dims(1) = 0
  do i = 1, hdf5_size
     dims(1) = dims(1) + sizes(i)
  enddo

  call h5screate_simple_f( rank, dims, filespace, ierr )
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_posr3p: h5screate_simple_f filespace dims:', dims
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
  call h5screate_simple_f( rank, sizes(hdf5_rank+1:), memspace, ierr )
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_posr3p: h5screate_simple_f memspace sizes:', sizes(hdf5_rank+1)
#endif

  ! 
  ! Select hyperslab in the file.
  !
  allocate(offset(hdf5_size))

  offset(1) = 0
  do i=2, hdf5_size
     offset(i) = offset(i-1) + sizes(i-1)
  enddo
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_posr3p: sizes', sizes, 'offset:', offset
#endif

  call h5dget_space_f( dset_id, filespace, ierr )
  call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, &
       offset(hdf5_rank+1:), sizes(hdf5_rank+1:), ierr )
#ifdef HDF5_DEBUG
  write(*,*) kfl_paral, 'hdf_posr3p: offset(hdf5_rank)', offset(hdf5_rank+1)
#endif

  !
  ! Write the dataset collectively. 
  !
#ifdef R4
  call h5dwrite_f( dset_id, H5T_NATIVE_REAL, buffer, sizes(hdf5_rank+1:), ierr, &
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#else
  call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, buffer, sizes(hdf5_rank+1:), ierr, &
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#endif


  deallocate(sizes)
  deallocate(buffer)
  deallocate(offset)

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

end subroutine hdf_posr3p
