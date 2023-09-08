!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     AUTHOR: Kenneth Leiter (kenneth.leiter@arl.army.mil)
!!
!!     Use the Xdmf Fortran Bindings to write out a simple mesh consisting of
!!     two hexahedrons.  Link against the XdmfFortran library to compile.
!!
!!     MODIFIED: Raúl de la Cruz (delacruz@bsc.es) April 2013
!!
!!     This example creates a H5 and XDMF files with a Temporal Collection
!!     Grid of a simple mesh composed of two hexahedrons with nodal and
!!     cell data.
!!
!!     HDF5 calls have been added in order to create an example to show
!!     how the new XdmfInterface class can be used in Fortran.
!!
!!     Please, check XdmfFortran.F90 interface file (Xdmf.mod) in Xdmf/libsrc/utils
!!     to review all the available XdmfInterface calls for Fortran.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM XdmfFortranExample
  USE XDMF
  USE HDF5 ! This module contains all necessary modules
  IMPLICIT NONE

  ! XDMF stuff
  INTEGER(8)                   :: obj
  INTEGER(8)                   :: ref_topo, ref_geom, ref_info, ref_attr1, ref_attr2
  CHARACTER(LEN=5)             :: group                    ! Full group name
  CHARACTER(LEN=9), PARAMETER  :: filename  = "my_output"  ! File name
  CHARACTER(LEN=8), PARAMETER  :: toponame  = "Topology"   ! Topology name
  CHARACTER(LEN=8), PARAMETER  :: geomname  = "Geometry"   ! Geometry name
  CHARACTER(LEN=4), PARAMETER  :: groupname = "Iter"       ! Group name
  CHARACTER(LEN=10), PARAMETER :: nodename  = "NodeValues" ! Node name
  CHARACTER(LEN=10), PARAMETER :: cellname  = "CellValues" ! Cell name

  ! Grid information
  INTEGER(4)     :: myConnections(8,2)
  REAL(4)        :: myPoints(3,3,4)
  REAL(8)        :: myCellAttribute(2), myNodeAttribute(3,4)

  ! HDF5 stuff
  INTEGER(HID_T) :: file_id       ! File identifier 
  INTEGER(HID_T) :: group_id      ! Group identifier 
  INTEGER(HID_T) :: dset_id       ! Dataset identifier 
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

  INTEGER(4)                     :: i, ierr, rank ! Error flags, Dataset rank 
  INTEGER(HSIZE_T), DIMENSION(3) :: dimsf         ! Dataset dimensions.


  ! Grid geometry
  myPoints(1,1,1) = 0
  myPoints(2,1,1) = 0
  myPoints(3,1,1) = 1
  myPoints(1,2,1) = 1
  myPoints(2,2,1) = 0
  myPoints(3,2,1) = 1
  myPoints(1,3,1) = 3
  myPoints(2,3,1) = 0
  myPoints(3,3,1) = 2
  myPoints(1,1,2) = 0
  myPoints(2,1,2) = 1
  myPoints(3,1,2) = 1
  myPoints(1,2,2) = 1
  myPoints(2,2,2) = 1
  myPoints(3,2,2) = 1
  myPoints(1,3,2) = 3
  myPoints(2,3,2) = 2
  myPoints(3,3,2) = 2
  myPoints(1,1,3) = 0
  myPoints(2,1,3) = 0
  myPoints(3,1,3) = -1
  myPoints(1,2,3) = 1
  myPoints(2,2,3) = 0
  myPoints(3,2,3) = -1
  myPoints(1,3,3) = 3
  myPoints(2,3,3) = 0
  myPoints(3,3,3) = -2
  myPoints(1,1,4) = 0
  myPoints(2,1,4) = 1
  myPoints(3,1,4) = -1
  myPoints(1,2,4) = 1
  myPoints(2,2,4) = 1
  myPoints(3,2,4) = -1
  myPoints(1,3,4) = 3
  myPoints(2,3,4) = 2
  myPoints(3,3,4) = -2

  ! Grid topology
  myConnections(1,1) = 0
  myConnections(2,1) = 1
  myConnections(3,1) = 7
  myConnections(4,1) = 6
  myConnections(5,1) = 3
  myConnections(6,1) = 4
  myConnections(7,1) = 10
  myConnections(8,1) = 9
  myConnections(1,2) = 1
  myConnections(2,2) = 2
  myConnections(3,2) = 8
  myConnections(4,2) = 7
  myConnections(5,2) = 4
  myConnections(6,2) = 5
  myConnections(7,2) = 11
  myConnections(8,2) = 10

  ! Node attributes
  myNodeAttribute(1,1) = 100
  myNodeAttribute(1,2) = 300
  myNodeAttribute(1,3) = 300
  myNodeAttribute(1,4) = 500
  myNodeAttribute(2,1) = 200
  myNodeAttribute(2,2) = 400
  myNodeAttribute(2,3) = 400
  myNodeAttribute(2,4) = 600
  myNodeAttribute(3,1) = 300
  myNodeAttribute(3,2) = 500
  myNodeAttribute(3,3) = 500
  myNodeAttribute(3,4) = 700

  ! Cell attributes
  myCellAttribute(1) = 100
  myCellAttribute(2) = 200

  ! Set references to null
  ref_topo  = 0
  ref_geom  = 0
  ref_info  = 0
  ref_attr1 = 0
  ref_attr2 = 0


  !
  ! Write HDF5 stuff
  !
  ! Initialize FORTRAN predefined datatypes
  CALL h5open_f(ierr)

  ! Create the file collectively.
  CALL h5fcreate_f(TRIM(filename)//'.h5', H5F_ACC_TRUNC_F, file_id, ierr)


  !
  ! Topology
  !

  ! Create the data space for the dataset.
  rank = 2
  dimsf(1) = 8
  dimsf(2) = 2
  CALL h5screate_simple_f(rank, dimsf, dspace_id, ierr)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(file_id, toponame, H5T_NATIVE_INTEGER, dspace_id, &
                   dset_id, ierr)

  ! Write the dataset collectively. 
  CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, myConnections, dimsf, ierr)

  ! Close dataspace.
  CALL h5sclose_f(dspace_id, ierr)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, ierr)


  !
  ! Geometry
  !

  ! Create the data space for the dataset.
  rank = 3
  dimsf(1) = 3
  dimsf(2) = 3
  dimsf(3) = 4
  CALL h5screate_simple_f(rank, dimsf, dspace_id, ierr)

  ! Create the dataset with default properties.
  CALL h5dcreate_f(file_id, geomname, H5T_NATIVE_REAL, dspace_id, &
                   dset_id, ierr)

  ! Write the dataset collectively. 
  CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, myPoints, dimsf, ierr)

  ! Close dataspace.
  CALL h5sclose_f(dspace_id, ierr)

  ! Close the dataset.
  CALL h5dclose_f(dset_id, ierr)


  do i= 1, 3
    ! Create group for node and cell data
    write(group, '(A,I1)') groupname, i
    CALL h5gcreate_f(file_id, group, group_id, ierr)


    !
    ! Node Attributes
    !

    ! Create the data space for the dataset.
    rank = 2
    dimsf(1) = 3
    dimsf(2) = 4
    CALL h5screate_simple_f(rank, dimsf, dspace_id, ierr)

    ! Create the dataset with default properties.
    CALL h5dcreate_f(group_id, nodename, H5T_NATIVE_DOUBLE, dspace_id, &
                     dset_id, ierr)

    ! Write the dataset collectively. 
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, myNodeAttribute, dimsf, ierr)

    ! Close dataspace.
    CALL h5sclose_f(dspace_id, ierr)

    ! Close the dataset.
    CALL h5dclose_f(dset_id, ierr)


    !
    ! Cell Attributes
    !

    ! Create the data space for the dataset.
    rank = 1
    dimsf(1) = 2
    CALL h5screate_simple_f(rank, dimsf, dspace_id, ierr)

    ! Create the dataset with default properties.
    CALL h5dcreate_f(group_id, cellname, H5T_NATIVE_DOUBLE, dspace_id, &
                     dset_id, ierr)

    ! Write the dataset collectively. 
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, myCellAttribute, dimsf, ierr)

    ! Close dataspace.
    CALL h5sclose_f(dspace_id, ierr)

    ! Close the dataset.
    CALL h5dclose_f(dset_id, ierr)

    ! Close the group.
    CALL h5gclose_f(group_id, ierr)

    myNodeAttribute = myNodeAttribute * 2
    myCellAttribute = myCellAttribute * 2
  enddo

  ! Close the file.
  CALL h5fclose_f(file_id, ierr)

  ! Close FORTRAN predefined datatypes.
  CALL h5close_f(ierr)


  !
  ! Write XDMF stuff
  !
  CALL XDMFINIT(obj, C_STR(filename), XDMF_TRUE)
  CALL XDMFSETDEBUG(obj, 1)

  ierr = XDMFADDGRID(obj, C_STR('Collection'), XDMF_GRID_COLLECTION, &
                     XDMF_GRID_COLLECTION_TEMPORAL)


  do i= 1, 3
    ! Grid
    ref_info  = 0
    ref_attr1 = 0
    ref_attr2 = 0
    write(group, '(A,I1)') groupname, i

    ierr = XDMFADDGRID(obj, C_STR('TestGrid'), XDMF_GRID_UNIFORM, XDMF_GRID_COLLECTION_UNSET)
    if (XDMFADDINFORMATION(obj, ref_info, C_STR('Info'), C_STR('This is a test')) == XDMF_FAIL) then
      STOP 'Error Adding Information Tag'
    endif
    if (XDMFSETGRIDTOPOLOGY(obj, ref_topo, XDMF_HEX, XDMF_INT32_TYPE, &
                            1, 2_8, 1, 16_8, C_STR(filename//'.h5:/'//toponame)) == XDMF_FAIL) then
      STOP 'Error Setting Grid Topology'
    endif
    if (XDMFSETGRIDGEOMETRY(obj, ref_geom, XDMF_GEOMETRY_XYZ, XDMF_FLOAT32_TYPE, &
                            1, 12_8, C_STR(filename//'.h5:/'//geomname)) == XDMF_FAIL) then
      STOP 'Error Setting Grid Geometry'
    endif

    CALL XDMFSETTIME(obj, REAL(i,8))
    if (XDMFADDGRIDATTRIBUTE(obj, ref_attr1, C_STR('NodeValues'),XDMF_FLOAT64_TYPE, XDMF_ATTRIBUTE_CENTER_NODE, &
                             XDMF_ATTRIBUTE_TYPE_SCALAR, 1, 12, C_STR(filename//'.h5:/'//group//'/'//nodename)) == XDMF_FAIL) then
      STOP 'Error Setting Grid Attribute'
    endif
    if (XDMFADDGRIDATTRIBUTE(obj, ref_attr2, C_STR('CellValues'), XDMF_FLOAT64_TYPE, XDMF_ATTRIBUTE_CENTER_CELL, &
                             XDMF_ATTRIBUTE_TYPE_SCALAR, 1, 2, C_STR(filename//'.h5:/'//group//'/'//cellname)) == XDMF_FAIL) then
      STOP 'Error Setting Grid Attribute'
    endif
    if (XDMFWRITEGRID(obj) == XDMF_FAIL) then
      STOP 'Error Building Grid Xml Structure'
    endif
    if (XDMFCLOSEGRID(obj) == XDMF_FAIL) then
      STOP 'Error Closing Grid'
    endif
  enddo

  ! Close Temporal Collection Grid
  ierr = XDMFCLOSEGRID(obj)

  ! Write output
  CALL XDMFWRITETOFILE(obj)
  !CALL XDMFSERIALIZE(obj)
  CALL XDMFCLOSE(obj)
END PROGRAM
