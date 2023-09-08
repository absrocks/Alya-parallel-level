module def_hdfpos

  !-----------------------------------------------------------------------
  !****f* hdfpos/def_hdfpos
  ! NAME
  !    def_hdfpos
  ! DESCRIPTION
  !    Heading for the HDF5 service
  ! USED BY
  !    Almost all hdf5 subroutines
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_domain
  use iso_c_binding
  use hdf5

  !------------------------------------------------------------------------
  ! Types
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! File names
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Initialization
  !------------------------------------------------------------------------

  ! Comunication variables
  integer(4)                :: &
       hdf5_group,             &   ! Group for HDF5 tasks
       hdf5_comm,              &   ! MPI Communicator for HDF5 tasks
       hdf5_rank,              &   ! Rank for HDF5 communicator
       hdf5_size,              &   ! Size for HDF5 communicator
       hdf5_info,              &   ! MPI Info identificator for HDF5 (hints)
       hdf5_ltype2xdmf(nelty), &   ! Array to translate from alya type to xdmf
       hdf5_ltype2vtk(nelty),  &   ! Array to translate from alya type to vtk
       hdf5_conne2vtk(nelty)       ! Array to translate from alya type to vtk

  integer(hsize_t)          :: &
       hdf5_connedim,          &    ! Connectivity HDF5 array dimension (XDMF Type+Nodes)
       hdf5_condifin,          &    ! Final connectivity sum for master
       hdf5_i                       ! Hdfd flag

  integer(hsize_t), pointer :: &
       hdf5_npoins(:),         &   ! npoin offset per task
       hdf5_nelems(:),         &   ! nelem offset per task
       hdf5_nvort(:),          &   ! number of vortex offset per task
       hdf5_filter(:)              ! number of points filtered per task

  logical                   :: &
       hdf5_opened,            &   ! Flag that indicates if opened file
       hdf5_mesh,              &   ! Flag that indicates if mesh dsm file
       hdf5_mesh_filter,       &   ! Flag that indicates if filter
       hdf5_group_filter,      &   ! Flag that indicates group of filter
       hdf5_restar_write,      &   ! Flag that indicates restart write
       hdf5_restar_read            ! Flag that indicates restart read

  ! Variables to manage h5 data files for each timestep
  integer(hid_t)            :: &
       hdf5_dcplist,           &    ! Data creation property list
       hdf5_dxplist,           &    ! Data transfer property list
       hdf5_faplist,           &    ! File access property list
       hdf5_fileid,            &    ! Collective file id
       group_id,               &    ! Group identifier
       hdf5_fileid_restart          ! Collective file id for restart

  character(150)            :: &
       hdf5_fil_error               ! File name for HDF5 output error

  integer(c_intptr_t)       :: &
       xdmf_obj                    ! C pointer to Xdmf object


end module def_hdfpos
