subroutine hdf_initia()
  use def_kintyp
  use def_hdfpos
  use def_parall
  use def_master
  use def_elmtyp
  use mod_parall, only : PAR_INTEGER,PAR_COMM_MY_CODE
  implicit none
#ifdef MPI_OFF
#else
  include 'mpif.h'
#endif

  integer(4)      :: ierr
  character        :: dummy
  integer(4)       :: stat, stats(13)
  integer(ip)      :: world_group, exclude
  integer(hsize_t) :: threshold, bsize
  character(len=10) :: bsizestr

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  ! Initialize open file flag
  hdf5_opened       = .false. 
  hdf5_mesh         = .false.
  hdf5_mesh_filter  = .false.
  hdf5_group_filter = .false.
  hdf5_restar_write = .false.
  hdf5_restar_read  = .false.

  ! Initialize hdf5_i : indice of time step
  hdf5_i = -1.0_ip


#ifdef MPI_OFF
  hdf5_rank = 0
  hdf5_size = 1

  !
  ! Initialize HDF5 FORTRAN predefined datatypes
  !
  !call h5open_f( ierr )

  ! Set access properties by default
  hdf5_dcplist = H5P_DEFAULT_F
  hdf5_dxplist = H5P_DEFAULT_F
  hdf5_faplist = H5P_DEFAULT_F
#else
  !
  ! Create MPI-IO communicator (Master not included)
  !

  ! Extract the original group handle
  call MPI_COMM_GROUP( PAR_COMM_MY_CODE, world_group, ierr )

  ! Remove master from the group
  exclude = 0
  call MPI_GROUP_EXCL( world_group, 1, exclude, hdf5_group, ierr )

  ! Create new communicator
  call MPI_COMM_CREATE( PAR_COMM_MY_CODE, hdf5_group, hdf5_comm, ierr)
  !hdf5_comm = PAR_COMM_MY_CODE

  ! Get my new rank/size for HDF5 communicator
  if (hdf5_comm .ne. MPI_COMM_NULL) then
    call MPI_COMM_RANK( hdf5_comm, hdf5_rank, ierr )
    call MPI_COMM_SIZE( hdf5_comm, hdf5_size, ierr )

#ifdef HDF5_DEBUG
    write(*,*) kfl_paral, 'hdf5_comm:', hdf5_comm, PAR_COMM_MY_CODE
    write(*,*) kfl_paral, 'hdf5_rank:', hdf5_rank
    write(*,*) kfl_paral, 'hdf5_size:', hdf5_size
#endif

    !
    ! Initialize HDF5 FORTRAN predefined datatypes
    !
    call h5open_f( ierr )
    if (ierr .lt. 0) then
      write(*,*) 'ERROR: Cannot initialize HDF5 library'
      return
    endif

 
! Default block size definition (MN3 GPFS: scratch and projects)
#define BLOCK_SIZE (4 * 1024 * 1024)

    ! Obtain preferred block size for the current File System
    ierr = STAT( adjustl(trim(namda))//'.dat', stats )
    if (ierr .eq. 0) then
      bsize = stats(12)
    else
      bsize = BLOCK_SIZE
    endif
    write (bsizestr,'(I10)') bsize

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! MPI-IO/GPFS Tuning: This may be hardware-dependent
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Improve Collective Read/Write through MPI-IO Hints (on research)
    ! On GPFS file systems IBM_io_buffer_size, IBM_largeblock_io and
    ! IBM_sparse_access might improve the performance for large
    ! number of MPI tasks and HDF5 module enabled with original write pattern
    call MPI_INFO_CREATE( hdf5_info, ierr )
    if (ierr .ne. MPI_SUCCESS) then
       write(*,*) 'Error creating MPI info'
    endif
    call MPI_INFO_SET( hdf5_info, 'IBM_io_buffer_size', '4M', ierr )
    if (ierr .ne. MPI_SUCCESS) then
      write(*,*) 'Error setting MPI-IO info hint IBM_io_buffer_size'
    endif
    call MPI_INFO_SET( hdf5_info, 'IBM_largeblock_io', 'true', ierr )
    if (ierr .ne. MPI_SUCCESS) then
      write(*,*) 'Error setting MPI-IO info hint IBM_largeblock_io'
    endif
    call MPI_INFO_SET( hdf5_info, 'IBM_sparse_access', 'false', ierr )
    if (ierr .ne. MPI_SUCCESS) then
      write(*,*) 'Error setting MPI-IO info hint IBM_sparse_access'
    endif
    call MPI_INFO_SET( hdf5_info, "striping_unit",bsizestr, ierr )
    if (ierr .ne. MPI_SUCCESS) then
      write(*,*) 'Error setting MPI-IO info hint striping_unit'
    endif
    call MPI_INFO_SET( hdf5_info, "romio_cb_write","enable", ierr )
    if (ierr .ne. MPI_SUCCESS) then
      write(*,*) 'Error setting MPI-IO info hint striping_unit'
    endif
    
    !ESTE AFECTA, ESTUDIAR SU IMPACTO
    !call MPI_INFO_SET( hdf5_info, "striping_factor","64", ierr )
    !if (ierr .ne. MPI_SUCCESS) then
    !  write(*,*) 'Error setting MPI-IO info hint striping_unit'
    !endif
    
    !call MPI_INFO_SET( hdf5_info, "cb_nodes","4", ierr )
    !if (ierr .ne. MPI_SUCCESS) then
    !  write(*,*) 'Error setting MPI-IO info hint striping_unit'
    !endif
    !call MPI_INFO_SET( hdf5_info, "romio_no_indep_rw","true", ierr )
    !if (ierr .ne. MPI_SUCCESS) then
    !  write(*,*) 'Error setting MPI-IO info hint striping_unit'
    !endif
    !
    ! Setup file access property list with parallel I/O access.
    !
    call h5pcreate_f( H5P_FILE_ACCESS_F, hdf5_faplist, ierr )
    !     previous call with normal I/O mpi comm (H5fd)
    call h5pset_fapl_mpio_f( hdf5_faplist, hdf5_comm, hdf5_info, ierr )
    !     new call with improve I/O mpi comm (H5fddsm)
    !     call h5pset_fapl_dsm_f( hdf5_faplist, hdf5_comm, ierr )

    !
    ! Create property list for collective dataset write
    !
    call h5pcreate_f( H5P_DATASET_XFER_F, hdf5_dxplist, ierr )
    call h5pset_dxpl_mpio_f( hdf5_dxplist, H5FD_MPIO_COLLECTIVE_F, ierr )


    ! Set dataset alignment for file access
    threshold = 1
    call h5pset_alignment_f( hdf5_faplist, threshold, bsize, ierr )
    if (ierr .lt. 0) then
      write(*,*) 'Error setting alignment for HDF5 library', &
                 threshold, bsize, ierr
    endif

    ! Sets data transfer size buffer
    call h5pset_buffer_f( hdf5_dxplist, bsize, ierr )
    if (ierr .lt. 0) then
      write(*,*) 'Error setting the data transfer size to ', &
                 bsize
    endif

    ! Sets data sieve size buffer
    !call h5pset_sieve_buf_size_f( hdf5_faplist, bsize, ierr )
    !if (ierr .lt. 0) then
    !  write(*,*) 'Error setting the data sieve buffer size to ', &
    !             bsize, ierr
    !endif

    ! Set data layout as compact
    ! Store raw data in the dataset object header in file. This should only be used
    ! for datasets with small amounts of raw data. The raw data size limit is 64K
    ! (65520 bytes). Attempting to create a dataset with raw data larger than this
    ! limit will cause the H5Dcreate call to fail. 
    !QUIZAS ESTA LINEA SOBRA, PORQUE SOBREESCRIBE LOS DEFAULTS O NO?
    call h5pcreate_f( H5P_DATASET_CREATE_F, hdf5_dcplist, ierr )
   !call h5pset_layout_f( hdf5_dcplist, H5D_COMPACT_F, ierr )
   !if (ierr .lt. 0) then
   !  write(*,*) 'Error setting the data layout to H5D_COMPACT_F ', ierr
   !endif

  ! Sanity check
  else if (kfl_paral .ne. 0) then
    call runend('Error assigning communicator for HDF5 service')
  endif
#endif


  !
  ! Allocate HDF5 arrays
  !
  allocate(hdf5_npoins(hdf5_size))
  allocate(hdf5_nelems(hdf5_size))


  ! HDF5 output error file name
  hdf5_fil_error = 'hdfpos.err'


  ! Conversion from Alya elmtype to XDMF standard
  hdf5_ltype2xdmf(POINT) = z'01'
  hdf5_ltype2xdmf(BAR02) = z'02' ! Implemented as a group of line segments
  hdf5_ltype2xdmf(BAR03) = z'02' ! Implemented as a group of line segments
  hdf5_ltype2xdmf(BAR04) = z'02' ! Implemented as a group of line segments
  hdf5_ltype2xdmf(TRI03) = z'04'
  hdf5_ltype2xdmf(TRI06) = z'24'
  hdf5_ltype2xdmf(QUA04) = z'05'
  hdf5_ltype2xdmf(QUA08) = z'25'
  hdf5_ltype2xdmf(QUA09) = z'25' ! Use QUA08 (does not exist in XDMF)
  hdf5_ltype2xdmf(QUA16) = z'25' ! Use QUA08 (does not exist in XDMF)
  hdf5_ltype2xdmf(TET04) = z'06'
  hdf5_ltype2xdmf(TET10) = z'26'
  hdf5_ltype2xdmf(PYR05) = z'07'
  hdf5_ltype2xdmf(PYR14) = z'27' ! Use PYR13 (does not exist in XDMF)
  hdf5_ltype2xdmf(PEN06) = z'08'
  hdf5_ltype2xdmf(PEN15) = z'28'
  hdf5_ltype2xdmf(PEN18) = z'29'
  hdf5_ltype2xdmf(HEX08) = z'09'
  hdf5_ltype2xdmf(HEX20) = z'30'
  hdf5_ltype2xdmf(HEX27) = z'32'
  hdf5_ltype2xdmf(HEX64) = z'32' ! Use HEX27 (does not exist in XDMF)

  ! Conversion from Alya elmtype to VTK standard
  
  hdf5_conne2vtk(TET04) = '4'
  hdf5_conne2vtk(PYR05) = '5'
  hdf5_conne2vtk(PEN06) = '6'
  hdf5_conne2vtk(HEX08) = '8'

  hdf5_ltype2vtk(TET04) = '10' ! VTK_TETRA
  hdf5_ltype2vtk(PYR05) = '14' ! VTK_PYRAMID
  hdf5_ltype2vtk(PEN06) = '13' ! VTK_WEDGE
  hdf5_ltype2vtk(HEX08) = '12' ! VTK_HEXAHEDRON

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine hdf_initia

