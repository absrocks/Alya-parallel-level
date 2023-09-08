subroutine hdf_posvor(itask)
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
  integer(4)                :: i, j, ierr, rank = 2,tot
  integer(hid_t)            :: filespace
  integer(hid_t)            :: dset_id     ! Dataset identifier
  integer(hid_t)            :: memspace    ! Memspace identifier
  integer(hsize_t), pointer :: offset(:,:)
  integer(hsize_t)          :: dims(2)
  integer(ip),      pointer :: sizee(:)
  
#ifdef EVENT
  call mpitrace_user_function(1)
#endif


  select case ( itask )

  case ( 0_ip )

     !
     ! Gather nvort offsets
     !
     allocate(sizee(hdf5_size))
     allocate(hdf5_nvort(hdf5_size))
     nvort_total=0_ip

#ifdef MPI_OFF
#else
     call MPI_ALLGATHER( nvort, 1, PAR_INTEGER, &
          sizee, 1, PAR_INTEGER, hdf5_comm, ierr )
#endif

     hdf5_nvort(1) = 0
     do i=2, hdf5_size
        hdf5_nvort(i) = hdf5_nvort(i-1) + sizee(i-1)
     enddo
     do i=1,hdf5_size
        nvort_total = nvort_total + sizee(i)
     end do

#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_posvor: hdf5_nvort', hdf5_nvort
     write(*,*) kfl_paral, 'hdf_posvor: nvort_total', nvort_total
#endif


     deallocate(sizee)


  case ( 1_ip )   
     !
     ! Create the data space for the dataset. 
     !
     dims(1) = ndime
     dims(2) = nvort_total
     !
     ! case of nvort=0
     !     
     if (dims(2)==0) then
        call h5screate_f(H5S_NULL_F, filespace,ierr)
     else
        call h5screate_simple_f( rank, dims, filespace, ierr )
     end if
     

#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_posvor: h5screate_simple_f filespace dims:', dims
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
     dims(1) = ndime
     dims(2) = nvort
     !
     ! case of nvort=0
     !     
     if (dims(2)==0) then
        call h5screate_f(H5S_NULL_F, memspace,ierr)
     else
        call h5screate_simple_f( rank, dims, memspace, ierr )
     end if

#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_posvor: h5screate_simple_f memspace dims:', dims
#endif

     ! 
     ! Select hyperslab in the file.
     !
     tot=0
     allocate(offset(2,hdf5_size))
     do i= 1, hdf5_size
        offset(1,i) = 0
        offset(2,i) = hdf5_nvort(i)
        tot=tot+hdf5_nvort(i)
     enddo
     call h5dget_space_f( dset_id, filespace, ierr )
     !
     ! case of nvort=0
     !     
     if (tot /= 0) then
        call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, offset(1:2,hdf5_rank+1), dims, ierr )
     end if
     
     
#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_posvor: offset(hdf5_rank+1)', offset(1:2,hdf5_rank+1)
     write(*,*) "hdf5_size",hdf5_size
     write(*,*) "hdf5_vorti",hdf5_nvort
     write(*,*) "kfl_paral, nvort",kfl_paral,nvort
#endif
     deallocate(offset)
     deallocate(hdf5_nvort)

     !
     ! Write the dataset collectively. 
     !
#ifdef R4
     call h5dwrite_f( dset_id, H5T_NATIVE_REAL, gevec_hdf, dims, ierr, &
          file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#else
     call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, gevec_hdf, dims, ierr, &
          file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#endif
#ifdef HDF5_DEBUG
     do i=1, nvort
        write(*,*) kfl_paral, 'hdf_posvor gevec_hdf:', gevec_hdf(1:SIZE(gevec_hdf, DIM=1), i)
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


  case ( 2_ip )
     !
     ! MPI_REDUCE : add nvort for each slave and put in nvort_total in master
     !

#ifdef MPI_OFF
#else
     call MPI_REDUCE( nvort, nvort_total, 1, PAR_INTEGER, &
                      MPI_SUM, 0, PAR_COMM_MY_CODE, ierr )
#endif

#ifdef HDF5_DEBUG
     write(*,*) kfl_paral, 'hdf_posvor: nvort_total', nvort_total
#endif


  end select


#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine hdf_posvor
