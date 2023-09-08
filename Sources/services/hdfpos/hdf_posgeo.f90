subroutine hdf_posgeo(itask)
  use def_kintyp
  use def_master
  use def_hdfpos
  use def_parall
  use def_domain
  use mod_parall, only : PAR_INTEGER,PAR_COMM_MY_CODE
!
!  use h5fddsm
!
  implicit none
#ifdef MPI_OFF
#else
  include 'mpif.h'
#endif

  integer(ip), intent(in)   :: itask

  integer(4)                :: i, ielem, inode, ibuf, ierr, rank1 = 1, rank2 = 2
  integer(hid_t)            :: file_id, filespace
  integer(hid_t)            :: dset_id     ! Dataset identifier
  integer(hid_t)            :: memspace    ! Memspace identifier
  integer(hsize_t)          :: dims(2)
  integer(hsize_t), pointer :: offset(:,:)
  integer(hsize_t), pointer :: sizes(:)
  integer(ip),       pointer :: buffer(:)
  integer(ip),      pointer :: sizee(:)
  
  TYPE(C_PTR) :: f_ptr
  INTEGER(HID_T) :: h5_kind_type_r, h5_kind_type_i ! HDF type corresponding to the specified KIND
  

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  select case ( itask )

  case ( 0_ip )
    !
    ! Gather npoin and nelem offsets
    !
    allocate(sizee(hdf5_size))
#ifdef MPI_OFF
#else
    call MPI_ALLGATHER( npoin, 1, PAR_INTEGER, &
                        sizee, 1, PAR_INTEGER, hdf5_comm, ierr )
#endif

    hdf5_npoins(1) = 0
    do i=2, hdf5_size
      hdf5_npoins(i) = hdf5_npoins(i-1) + sizee(i-1)
    enddo

#ifdef MPI_OFF
#else
    call MPI_ALLGATHER( nelem, 1, PAR_INTEGER, &
                        sizee, 1, PAR_INTEGER, hdf5_comm, ierr )
#endif

    hdf5_nelems(1) = 0
    do i=2, hdf5_size
      hdf5_nelems(i) = hdf5_nelems(i-1) + sizee(i-1)
    enddo

#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: hdf5_npoins', hdf5_npoins
    write(*,*) kfl_paral, 'hdf_posgeo: hdf5_nelems', hdf5_nelems
#endif
    deallocate(sizee)

  case ( 1_ip )
    !
    ! POSTPROCESS COORDINATES
    !
    wopos_hdf(1) = 'COORD'

    !
    ! Create the file collectively.
    !
    call h5fcreate_f( fil_postp, H5F_ACC_TRUNC_F, file_id, ierr, access_prp=hdf5_faplist )

    !
    ! Create the data space for the dataset. 
    !
    dims(1) = ndime
    dims(2) = npoin_total
    call h5screate_simple_f( rank2, dims, filespace, ierr )
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: h5screate_simple_f filespace dims:', dims
#endif

    !
    ! Create the dataset with default properties.
    !
#ifdef R4
    call h5dcreate_f( file_id, wopos_hdf(1), H5T_NATIVE_REAL, filespace, &
                      dset_id, ierr, dcpl_id = hdf5_dcplist )
#else
    call h5dcreate_f( file_id, wopos_hdf(1), H5T_NATIVE_DOUBLE, filespace, &
                      dset_id, ierr, dcpl_id = hdf5_dcplist )
#endif
    call h5sclose_f( filespace, ierr )

    !
    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file. 
    !
    dims(1) = ndime
    dims(2) = npoin
    call h5screate_simple_f( rank2, dims, memspace, ierr )
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: h5screate_simple_f memspace dims:', dims
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
    write(*,*) kfl_paral, 'hdf_posgeo: offset(hdf5_rank+1)', offset(1:2,hdf5_rank+1)
#endif
    deallocate(offset)

    !
    ! Write the dataset collectively. 
    !
#ifdef R4
    call h5dwrite_f( dset_id, H5T_NATIVE_REAL, coord, dims, ierr, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#else
    call h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, coord, dims, ierr, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = hdf5_dxplist )
#endif
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: coord:', SIZE(coord), ' dims:', dims
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


    !
    ! POSTPROCESS CONNECTIVITIES
    !
    wopos_hdf(1) = 'CONNE'

    !
    ! Create the data space for the dataset. 
    !

    ! Calculate the required size for the buffer
    !
    ! Mixed topology is used in HDF5/XDMF format
    ! Needs space for 'element type + element nodes'

    allocate(sizes(hdf5_size))
    ! Add one entry per element to store element type
    hdf5_connedim = nelem
    do ielem = 1, nelem
      hdf5_connedim = hdf5_connedim + nnode(ltype(ielem))
    enddo
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: hdf5_connedim:', hdf5_connedim
#endif

    ! Copy XDMF type + lnods data in a sequential buffer
    ! XDMF-ElementType + Nodes (in 0-BaseOffset)
    allocate(buffer(hdf5_connedim))
    ibuf = 1
    do ielem = 1, nelem
      buffer(ibuf) = hdf5_ltype2xdmf(ltype(ielem))
      do inode = 1, nnode(ltype(ielem))
        buffer(ibuf + inode) = lnods(inode,ielem) + hdf5_npoins(hdf5_rank+1) - 1
      enddo
      ibuf = ibuf + 1 + nnode(ltype(ielem))
    enddo
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: buffer:', buffer, 'ibuf:', ibuf
#endif

    ! Sanity check for memory copy
    if ((ibuf-1) .ne. hdf5_connedim) then
      call runend('Error in hdf_posgeo copying lnods structure to sequential buffer')
    endif

#ifdef MPI_OFF
#else
  call MPI_ALLGATHER( hdf5_connedim, 1, MPI_INTEGER8, &
                      sizes, 1, MPI_INTEGER8, hdf5_comm, ierr )
#endif

    dims(1) = 0
    do i = 1, hdf5_size
      dims(1) = dims(1) + sizes(i)
    enddo

    call h5screate_simple_f( rank1, dims, filespace, ierr )
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: h5screate_simple_f filespace dims:', dims
#endif

    !
    ! Create the dataset with default properties.
    !
#ifdef I8
    ! SOLVED, WORKS WITH FORTRAN 2003 AND INT 8
    ! This should be H5T_NATIVE_LLONG, however this type is not supported in
    ! Fortran interface, therefore H5T_NATIVE_INTEGER is used (Loss of ! precision).
    ! Check: New Features in the HDF5 Fortran Library: Adding support for the Fortran 2003 Standard
    !        http://comments.gmane.org/gmane.comp.programming.hdf/3467
    !write(*,*) 'WARNING: h5dcreate_f is using H5T_NATIVE_INTEGER instead of H5T_NATIVE_LLONG'
    h5_kind_type_i = h5kind_to_type(ip,H5_INTEGER_KIND)    
    call h5dcreate_f( file_id, wopos_hdf(1), h5_kind_type_i, filespace, &
                      dset_id, ierr, dcpl_id = hdf5_dcplist )
#else
    call h5dcreate_f( file_id, wopos_hdf(1), H5T_NATIVE_INTEGER, filespace, &
                      dset_id, ierr, dcpl_id = hdf5_dcplist )
#endif
    call h5sclose_f( filespace, ierr )

    !
    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file. 
    !
    call h5screate_simple_f( rank1, sizes(hdf5_rank+1:), memspace, ierr )
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: h5screate_simple_f memspace sizes:', sizes(hdf5_rank+1)
#endif

    ! 
    ! Select hyperslab in the file.
    !
    allocate(offset(1,hdf5_size))

    offset(1,1) = 0
    do i=2, hdf5_size
      offset(1,i) = offset(1,i-1) + sizes(i-1)
    enddo
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: sizes', sizes, 'offset:', offset
#endif

    call h5dget_space_f( dset_id, filespace, ierr )
    call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, &
                                offset(1,hdf5_rank+1:), sizes(hdf5_rank+1:), ierr )
#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf_posgeo: offset(hdf5_rank)', offset(1,hdf5_rank+1)
#endif

    !
    ! Write the dataset collectively. 
    !
#ifdef I8
    ! This should be H5T_NATIVE_LLONG, however this type is not supported in
    ! Fortran interface, therefore H5T_NATIVE_INTEGER is used (Loss of ! precision).
    !write(*,*) 'WARNING: h5dwrite_f is using H5T_NATIVE_INTEGER instead of H5T_NATIVE_LLONG'
    f_ptr = C_LOC(buffer(1))
    CALL h5dwrite_f( dset_id, h5_kind_type_i, f_ptr, ierr, file_space_id = filespace, &
                     mem_space_id = memspace, xfer_prp = hdf5_dxplist)
#else
    call h5dwrite_f( dset_id, H5T_NATIVE_INTEGER, buffer, sizes(hdf5_rank+1:), ierr, &
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

    !
    ! Close the file.
    !
    call h5fclose_f( file_id, ierr )

  end select

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine hdf_posgeo
