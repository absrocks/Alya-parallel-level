subroutine hdf_restar()

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

  integer(4)                :: i, j, ierr, ibuf, rank1 = 1, rank2 = 2, ielem
  integer(hid_t)            :: filespace
  integer(hid_t)            :: dset_id     ! Dataset identifier
  integer(hid_t)            :: memspace    ! Memspace identifier
  integer(hsize_t)          :: dims1(1),dims2(2)
  integer(hsize_t), pointer :: offset(:,:)

  integer(hsize_t), pointer :: offset1(:)
  integer(hsize_t), pointer :: sizes(:)
  real(rp),         pointer :: buffer(:)
  integer(ip)               :: pgaus,pelty,igaus,idime,iprev


  if ( kfl_reawr == 0 ) return
  !
  ! do nothing ! 
  !
  if ( kfl_reawr == 1 ) then
     !
     ! Read the restart file
     !
     if (ISLAVE) then
        !
        ! call restart name
        !
        call hdf_openfi(5_ip)
        !
        ! Open an existing file.
        !
        call h5fopen_f (fil_rstar, H5F_ACC_RDONLY_F, hdf5_fileid, ierr)
        !
        ! Open an existing dataset. 
        !
        call h5dopen_f (hdf5_fileid, wopos_hdf(1), dset_id, ierr)
        !
        ! Read the dataset.
        !
        if      (wopos_hdf(2) == 'SCALA') then
           dims1(1) = npoin
           !
           ! Get dataset's dataspace identifier.
           !
           call h5dget_space_f( dset_id, filespace, ierr )          
           !
           ! Select hyperslab in the dataset.
           !
           call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, hdf5_npoins(hdf5_rank+1:), dims1, ierr )
           !
           ! Create memory dataspace.
           !
           call h5screate_simple_f( rank1, dims1, memspace, ierr )
           !
           ! Read data from hyperslab in the file into the hyperslab in memory 
           !
#ifdef R4
           call h5dread_f(dset_id, H5T_NATIVE_REAL, gesca_hdf, dims1, ierr, memspace, filespace )
#else
           call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, gesca_hdf, dims1, ierr, memspace, filespace )
#endif
           !
           !
           !
        else if (wopos_hdf(2) == 'VECTO') then
           dims2(1) = SIZE(gevec_hdf, DIM=1)
           dims2(2) = npoin
           !
           ! Get dataset's dataspace identifier.
           !
           call h5dget_space_f( dset_id, filespace, ierr )
           !
           ! Select hyperslab in the dataset.
           !
           allocate(offset(2,hdf5_size))
           do i= 1, hdf5_size
              offset(1,i) = 0
              offset(2,i) = hdf5_npoins(i)
           enddo
           call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, offset(1:2,hdf5_rank+1), dims2, ierr )
           deallocate(offset)
           !
           ! Create memory dataspace.
           !
           call h5screate_simple_f( rank2, dims2, memspace, ierr )
           !
           ! Read the data set
           !
#ifdef R4
           call h5dread_f(dset_id, H5T_NATIVE_REAL, gevec_hdf, dims2, ierr, memspace, filespace )
#else
           call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, gevec_hdf, dims2, ierr, memspace, filespace )
#endif
           !
           !
           ! 
        else if (wopos_hdf(2) == 'R3PVE') then
           !
           ! Calculate the required size for the buffer(only one time step dim1*dim2)
           !
           allocate(sizes(hdf5_size))
           dims1(1) = 0
           do ielem = 1, nelem
              dims1(1) = dims1(1) + SIZE(ger3p_hdf(ielem)%a, DIM=1)*SIZE(ger3p_hdf(ielem)%a, DIM=2)
           enddo
#ifdef MPI_OFF
#else
           call MPI_ALLGATHER( dims1, 1, MPI_INTEGER8, &
                sizes, 1, MPI_INTEGER8, hdf5_comm, ierr )
#endif

           allocate(offset1(hdf5_size))
           offset1(1) = 0
           do i=2, hdf5_size
              offset1(i) = offset1(i-1) + sizes(i-1)
           enddo
           !
           ! Copy r3p data in a sequential buffer
           !
           allocate(buffer(dims1(1)))!  
           ibuf = 1
           do ielem = 1, nelem    
              buffer(ibuf:ibuf+ SIZE(ger3p_hdf(ielem)%a, DIM=1)*SIZE(ger3p_hdf(ielem)%a, DIM=2)-1 ) & 
                   = RESHAPE(ger3p_hdf(ielem)%a, (/SIZE(ger3p_hdf(ielem)%a, DIM=1)*SIZE(ger3p_hdf(ielem)%a, DIM=2)/))
              ibuf = ibuf + SIZE(ger3p_hdf(ielem)%a, DIM=1)*SIZE(ger3p_hdf(ielem)%a, DIM=2)          
           enddo
           !
           ! Sanity check for memory copy
           !
           if ((ibuf-1) .ne. dims1(1)) then
              call runend('Error in hdf_restart(hdf_posr3p) copying r3p structure to sequential buffer')
           endif
           !
           ! Get dataset's dataspace identifier.
           !
           call h5dget_space_f( dset_id, filespace, ierr )
           !
           ! Select hyperslab in the dataset.
           !
           call h5sselect_hyperslab_f( filespace, H5S_SELECT_SET_F, offset1(hdf5_rank+1:), dims1, ierr )
           !
           ! Create memory dataspace.
           !
           call h5screate_simple_f( rank1, dims1, memspace, ierr )
           !
           ! Read the data set
           !
#ifdef R4           
           call h5dread_f(dset_id, H5T_NATIVE_REAL, buffer, dims1, ierr, memspace, filespace )
#else
           call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, dims1, ierr, memspace, filespace )
#endif
           i=0
           do ielem = 1,nelem
              pelty = abs(ltype(ielem))
              pgaus = ngaus(pelty)            
                 do igaus = 1,pgaus
                    do idime = 1,ndime
                       i=i+1
                       ger3p_hdf(ielem)%a(idime,igaus,2) = buffer(i)
                       !vesgs(ielem)%a(idime,igaus,2) = vesgs(ielem)%a(idime,igaus,1)
                    end do
                 end do
           end do
           !
           ! dealloc
           !
           deallocate(sizes)
           deallocate(buffer)
           deallocate(offset1)

        end if
        !
        ! Close the dataset.
        !
        CALL h5dclose_f(dset_id, ierr)
        !
        ! Close the file.
        !
        CALL h5fclose_f(hdf5_fileid, ierr)
     endif

  else if (kfl_reawr == 2) then
     !
     ! write the restart file
     !
     if (.not. hdf5_restar_write) then
        !
        ! Open H5 timestep restart file
        !
        hdf5_restar_write = .true.
        if (ISLAVE) then
           call hdf_openfi(4_ip)
        endif
     endif
     if (ISLAVE) then
        if      (wopos_hdf(2) == 'SCALA') then
           call hdf_possca()
        else if (wopos_hdf(2) == 'VECTO') then
           call hdf_posvec()
        else if (wopos_hdf(2) == 'R3PVE') then
           call hdf_posr3p()   
        endif
     end if
     
     
  end if



end subroutine hdf_restar
