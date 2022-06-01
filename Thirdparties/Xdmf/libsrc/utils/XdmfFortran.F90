!>*******************************************************************
!>                              XDMF
!>                  eXtensible Data Model and Format
!>
!> @file XdmfFortran.F90
!> @details Fortran 90 Module for XdmfInterface class
!>
!> Id : $Id $
!> Date : $Date $
!> Version : $Revision $
!>
!> @author
!>    Raul de la Cruz
!>    delacruz@bsc.es
!>    Barcelona Supercomputing Center - CNS
!>    Barcelona, Spain
!> @date
!>    February-July 2013                                              
!>                                                               
!>    This software is distributed WITHOUT ANY WARRANTY; without 
!>    even the implied warranty of MERCHANTABILITY or FITNESS    
!>    FOR A PARTICULAR PURPOSE.  See the above copyright notice  
!>    for more information.                                      
!>
!>*******************************************************************


!> @var This macro specifies that type and shape matching rules related
!> to explicit interfaces are to be ignored
!> We add comments in the MACRO call to avoid expanding MACRO to nothing
#ifdef __INTEL_COMPILER
# define NO_ARG_CHECK( var ) \
      DEC$ ATTRIBUTES NO_ARG_CHECK :: var !< Intel compiler (default)
#elif  __CRAY
# define NO_ARG_CHECK( var ) \
      DIR$ IGNORE_TKR var                 !< Cray compiler
#elif __IBM__
# define NO_ARG_CHECK( var ) \
      IBM* IGNORE_TKR var                 !< IBM compiler
#elif defined(__GFORTRAN__) && \
      (__GNUC__ >= 4 && (__GNUC_MINOR__ >= 9))
# define NO_ARG_CHECK( var ) \
      GCC$ ATTRIBUTES NO_ARG_CHECK :: var !< GNU compiler (from 4.9)
#else
# define NO_ARG_CHECK( var ) \
      $PRAGMA IGNORE_TKR var              !< Other compilers
#endif


!> This a Fortran 90 module wrapper for the XdmfInterface C class.
!> The XdmfInterface eases the Xdmf library utilization by offering
!> to the user a set of functions that Read/Write Xmdf (.xmf) files.
!> Under Examples/Cxx directory of Xdmf library you will find
!> some example of how to use this module in Fortran.
module Xdmf
  use iso_c_binding
  implicit none

  ! Version
  INTEGER(C_INT),   PARAMETER :: XDMF_VERSION_MAJOR  = 2      !< Major Xdmf version
  INTEGER(C_INT),   PARAMETER :: XDMF_VERSION_MINOR  = 1      !< Minor Xdmf version
  CHARACTER(LEN=*), PARAMETER :: XDMF_VERSION_STRING = '2.1'  !< Xdmf version

  ! Success and fail codes of Xdmf library
  INTEGER(C_INT), PARAMETER :: XDMF_SUCCESS =  1  !< Success Xdmf code
  INTEGER(C_INT), PARAMETER :: XDMF_FAIL    = -1  !< Fail Xdmf code

  ! True and false codes of Xdmf library
  INTEGER(C_INT), PARAMETER :: XDMF_TRUE  = 1 !< True Xdmf code
  INTEGER(C_INT), PARAMETER :: XDMF_FALSE = 0 !< False Xdmf code

  INTEGER(C_INT), PARAMETER :: XDMF_MAX_DIMENSION     = 10
  INTEGER(C_INT), PARAMETER :: XDMF_MAX_STRING_LENGTH = 1024

  ! Number Types
  INTEGER(C_INT), PARAMETER :: XDMF_UNKNOWN_TYPE  = -1  !< Unkown type
  INTEGER(C_INT), PARAMETER :: XDMF_INT8_TYPE     =  1  !< INTEGER*1 type
  INTEGER(C_INT), PARAMETER :: XDMF_INT16_TYPE    =  6  !< INTEGER*2 type
  INTEGER(C_INT), PARAMETER :: XDMF_INT32_TYPE    =  2  !< INTEGER*4 type
  INTEGER(C_INT), PARAMETER :: XDMF_INT64_TYPE    =  3  !< INTEGER*8 type
  INTEGER(C_INT), PARAMETER :: XDMF_FLOAT32_TYPE  =  4  !< REAL*4 type
  INTEGER(C_INT), PARAMETER :: XDMF_FLOAT64_TYPE  =  5  !< REAL*8 type
  INTEGER(C_INT), PARAMETER :: XDMF_UINT8_TYPE    =  7  !< UNSIGNED*1 type
  INTEGER(C_INT), PARAMETER :: XDMF_UINT16_TYPE   =  8  !< UNSIGNED*2 type
  INTEGER(C_INT), PARAMETER :: XDMF_UINT32_TYPE   =  9  !< UNSIGNED*4 type
  INTEGER(C_INT), PARAMETER :: XDMF_UINT64_TYPE   =  10 !< UNSIGNED*8 type
  INTEGER(C_INT), PARAMETER :: XDMF_COMPOUND_TYPE =  Z'10'  !< Compound type

  ! General Uniform Organization
  INTEGER(C_INT), PARAMETER :: XDMF_STRUCTURED   = 0  !< Structured mesh
  INTEGER(C_INT), PARAMETER :: XDMF_UNSTRUCTURED = 1  !< Unstructured mesh

  ! Topologies
  INTEGER(C_INT), PARAMETER :: XDMF_NOTOPOLOGY   = Z'0'    !< No topology
  INTEGER(C_INT), PARAMETER :: XDMF_POLYVERTEX   = Z'1'    !< Group of unconnected points
  INTEGER(C_INT), PARAMETER :: XDMF_POLYLINE     = Z'2'    !< Group of line segments
  INTEGER(C_INT), PARAMETER :: XDMF_POLYGON      = Z'3'    !< Points define a polygon
  INTEGER(C_INT), PARAMETER :: XDMF_TRI          = Z'4'    !< Triangle
  INTEGER(C_INT), PARAMETER :: XDMF_QUAD         = Z'5'    !< Quadrilateral
  INTEGER(C_INT), PARAMETER :: XDMF_TET          = Z'6'    !< Tetrahedron
  INTEGER(C_INT), PARAMETER :: XDMF_PYRAMID      = Z'7'    !< Pyramid
  INTEGER(C_INT), PARAMETER :: XDMF_WEDGE        = Z'8'    !< Wedge
  INTEGER(C_INT), PARAMETER :: XDMF_HEX          = Z'9'    !< Hexahedron
  INTEGER(C_INT), PARAMETER :: XDMF_EDGE_3       = Z'0022' !< Edge (Quadratic with 3 points)
  INTEGER(C_INT), PARAMETER :: XDMF_TRI_6        = Z'0024' !< Triangle (6 points)
  INTEGER(C_INT), PARAMETER :: XDMF_QUAD_8       = Z'0025' !< Quadrilateral (8 points)
  INTEGER(C_INT), PARAMETER :: XDMF_TET_10       = Z'0026' !< Tetrahedron (10 points)
  INTEGER(C_INT), PARAMETER :: XDMF_PYRAMID_13   = Z'0027' !< Pyramid (13 points)
  INTEGER(C_INT), PARAMETER :: XDMF_WEDGE_15     = Z'0028' !< Wedge (15 points)
  INTEGER(C_INT), PARAMETER :: XDMF_WEDGE_18     = Z'0029' !< Wedge (18 points)
  INTEGER(C_INT), PARAMETER :: XDMF_HEX_20       = Z'0030' !< Hexahedron (20 points)
  INTEGER(C_INT), PARAMETER :: XDMF_HEX_24       = Z'0031' !< Hexahedron (24 points)
  INTEGER(C_INT), PARAMETER :: XDMF_HEX_27       = Z'0032' !< Hexahedron (27 points)
  INTEGER(C_INT), PARAMETER :: XDMF_MIXED        = Z'0070' !< Mixture of unstructured cells
  INTEGER(C_INT), PARAMETER :: XDMF_2DSMESH      = Z'0100' !< 2D Curvilinear structured mesh
  INTEGER(C_INT), PARAMETER :: XDMF_2DRECTMESH   = Z'0101' !< 2D Rectilinear structured mesh (axis are perpendicular)
  INTEGER(C_INT), PARAMETER :: XDMF_2DCORECTMESH = Z'0102' !< 2D CoRectilinear structured mesh (axis are perpendicular and spacing is constant)
  INTEGER(C_INT), PARAMETER :: XDMF_3DSMESH      = Z'1100' !< 3D Curvilinear structured mesh
  INTEGER(C_INT), PARAMETER :: XDMF_3DRECTMESH   = Z'1101' !< 3D Rectilinear structured mesh (axis are perpendicular)
  INTEGER(C_INT), PARAMETER :: XDMF_3DCORECTMESH = Z'1102' !< 3D CoRectilinear structured mesh (axis are perpendicular and spacing is constant)

  ! Geometries
  INTEGER(C_INT), PARAMETER :: XDMF_GEOMETRY_NONE          = 0 !< No geometry
  INTEGER(C_INT), PARAMETER :: XDMF_GEOMETRY_XYZ           = 1 !< Interlaced locations (3D)
  INTEGER(C_INT), PARAMETER :: XDMF_GEOMETRY_XY            = 2 !< Interlaced locations (2D)
  INTEGER(C_INT), PARAMETER :: XDMF_GEOMETRY_X_Y_Z         = 3 !< X, Y and Z in separated arrays (3D)
  INTEGER(C_INT), PARAMETER :: XDMF_GEOMETRY_X_Y           = 4 !< X and Y in separated arrays (2D)
  INTEGER(C_INT), PARAMETER :: XDMF_GEOMETRY_VXVYVZ        = 5 !< Three arrays, one for each axis of size NX, NY, NZ dimensions (for 3DRECTMESH)
  INTEGER(C_INT), PARAMETER :: XDMF_GEOMETRY_ORIGIN_DXDYDZ = 6 !< Six values: Ox, Oy, Oz (origin) + Dx, Dy, Dz (discretization on each axis) (for 3DCORECTMESH)
  INTEGER(C_INT), PARAMETER :: XDMF_GEOMETRY_VXVY          = 7 !< Two arrays, one for each axis of size NX, NY dimensions (for 2DRECTMESH)
  INTEGER(C_INT), PARAMETER :: XDMF_GEOMETRY_ORIGIN_DXDY   = 8 !< Four values: Ox, Oy (origin) + Dx, Dy (discretization on each axis) (for 2DCORECTMESH)

  ! Grid collections
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_UNIFORM    = Z'00000' !< Homogeneous single grid (i.e. a pile of triangles)
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_COLLECTION = Z'10000' !< Array of Uniform Grids all with the same Attributes Grid
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_TREE       = Z'20000' !< Hierarchical Grid
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_SUBSET     = Z'40000' !< Portion of another Grid
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_UNSET      = Z'0FFFF' !< Grid collection type not set

  INTEGER(C_INT), PARAMETER :: XDMF_GRID_MASK       = Z'F0000' !< Used to check if Grid is Uniform or a Collection (Type xor XDMF_GRID_MASK = XdmfTopology Type)

  INTEGER(C_INT), PARAMETER :: XDMF_GRID_SECTION_ALL         = Z'100000'
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_SECTION_DATA_ITEM   = Z'200000'
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_SECTION_MASK        = Z'F00000'

  ! Types of Grid collections
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_COLLECTION_TEMPORAL = Z'0001'  !< Grid collection specifies an iteration in time of a Grid
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_COLLECTION_SPATIAL  = Z'0002'  !< Grid collection of a set of different Grids in same iteration of time
  INTEGER(C_INT), PARAMETER :: XDMF_GRID_COLLECTION_UNSET    = Z'0FFFF' !< Grid collection type not specified

  ! Value Types
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_TYPE_NONE     = 0 !< No type of data item specified
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_TYPE_SCALAR   = 1 !< Scalar type (1 value)
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_TYPE_VECTOR   = 2 !< Vector type (array of values)
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_TYPE_TENSOR   = 3 !< Tensor type (9 values)
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_TYPE_MATRIX   = 4 !< An arbitrary NxM matrix
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_TYPE_TENSOR6  = 5 !< Symmetrical tensor type (6 values) 
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_TYPE_GLOBALID = 6

  ! Where Values are Assigned
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_CENTER_GRID = 0 !< Centered to Grid (i.e. Material type)
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_CENTER_CELL = 1 !< Centered to a Cell
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_CENTER_FACE = 2 !< Centered to a Face
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_CENTER_EDGE = 3 !< Centered to a Edge
  INTEGER(C_INT), PARAMETER :: XDMF_ATTRIBUTE_CENTER_NODE = 4 !< Centered to a Node


  interface

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Initialize a new Xdmf file.
    !> @details
    !> Initializes a XdmfInterface object for a given output filename (.xdm)
    !> Supports generating HDF externally (default) or through Xdmf library
    !>
    !> @param[out] object       Variable containing the XdmfInterface object ID
    !> @param[in]  outputName   Name of the Xdmf file (.xdm)
    !> @param[in]  externalHDF  Specifies whether Xdmf writes HDF data
    !>                          (XDMF_FALSE) or not (XDMF_TRUE)
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfInit( object, outputName, externalHDF ) bind(c, name='XdmfInit')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(out)   :: object
      character(KIND=C_CHAR), intent(in) :: outputName(*)
      integer(C_INT), intent(in), value  :: externalHDF
    end subroutine XdmfInit

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Set Time.
    !> @details
    !> Sets the current Grid time, useful for temporal collection.
    !>
    !> @param[in]  object Variable containing the XdmfInterface object ID
    !> @param[in]  t      Time in double precision format (real*8)
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfSetTime( object, t ) bind(c, name='XdmfSetTime')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
      real(C_DOUBLE), intent(in), value      :: t
    end subroutine XdmfSetTime

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Get Time.
    !> @details
    !> Gets the current Grid time, useful for temporal collection.
    !>
    !> @param[in]  object Variable containing the XdmfInterface object ID
    !> @return Time in double precision format (real*8)
    !---------------------------------------------------------------------------  
    function XdmfGetTime( object ) bind(c, name='XdmfGetTime')
      use iso_c_binding
      implicit none
      real(C_DOUBLE)                         :: XdmfGetTime
      integer(C_INTPTR_T), intent(in), value :: object
    end function XdmfGetTime

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Set External HDF flag.
    !> @details
    !> Changes the status of the external HDF flag. If externalHDF = XDMF_TRUE,
    !> the Xdmf library writes the HDF data, otherwise (XDMF_FALSE) not.
    !>
    !> @param[in]  object       Variable containing the XdmfInterface object ID
    !> @param[in]  externalHDF  Value for the externalHDF flag (XDMF_TRUE/XDMF_FALSE)
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfSetExternalHDF( object, externalHDF ) bind(c, name='XdmfSetExternalHDF')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INT), intent(in), value      :: externalHDF
    end subroutine XdmfSetExternalHDF

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Get External HDF flag.
    !> @details
    !> Get the status of the external HDF flag. If externalHDF = XDMF_TRUE,
    !> the Xdmf library writes the HDF data, otherwise (XDMF_FALSE) not.
    !>
    !> @param[in]  object       Variable containing the XdmfInterface object ID
    !> @return Value for the externalHDF flag (XDMF_TRUE/XDMF_FALSE)
    !---------------------------------------------------------------------------  
    function XdmfGetExternalHDF( object ) bind(c, name='XdmfGetExternalHDF')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfGetExternalHDF
      integer(C_INTPTR_T), intent(in), value :: object
    end function XdmfGetExternalHDF

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Set Debug flag.
    !> @details
    !> Changes the status of the debug flag. If value = XDMF_TRUE,
    !> the Xdmf library writes debug output, otherwise (XDMF_FALSE) not.
    !>
    !> @param[in]  object       Variable containing the XdmfInterface object ID
    !> @param[in]  externalHDF  Value for the debug flag (XDMF_TRUE/XDMF_FALSE)
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfSetDebug( object, value ) bind(c, name='XdmfSetDebug')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INT), intent(in), value      :: value
    end subroutine XdmfSetDebug

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Get Debug flag.
    !> @details
    !> Get the status of the debug flag. If value = XDMF_TRUE,
    !> the Xdmf library writes debug output, otherwise (XDMF_FALSE) not.
    !>
    !> @param[in]  object       Variable containing the XdmfInterface object ID      
    !> @return Value for the debug flag (XDMF_TRUE/XDMF_FALSE)
    !---------------------------------------------------------------------------  
    function XdmfGetDebug( object ) bind(c, name='XdmfGetDebug')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfGetDebug
      integer(C_INTPTR_T), intent(in), value :: object
    end function XdmfGetDebug

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Add a XdmfGrid.
    !> @details
    !> Add a regular Grid or Collection to the XdmfDOM.
    !> Collections can be 'Spatial' or 'Temporal' type.
    !> Nested collections are supported.
    !>
    !> GridType can be: \n
    !>    XDMF_GRID_UNIFORM \n
    !>    XDMF_GRID_COLLECTION \n
    !>    XDMF_GRID_TREE \n
    !>    XDMF_GRID_SUBSET \n
    !>    XDMF_GRID_UNSET \n
    !>
    !> CollectionType can be: \n
    !>    XDMF_GRID_COLLECTION_TEMPORAL \n
    !>    XDMF_GRID_COLLECTION_SPATIAL \n
    !>    XMDF_GRID_COLLECTION_UNSET \n
    !>
    !> @param[in]  object         Variable containing the XdmfInterface object ID
    !> @param[in]  gridName       Name of the Grid to add
    !> @param[in]  gridType       Type of Grid (see GridType description above)
    !> @param[in]  collectionType Type of Collection if Grid is a Collection
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfAddGrid( object, gridName, gridType, collectionType ) bind(c, name='XdmfAddGrid')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfAddGrid
      integer(C_INTPTR_T), intent(in), value :: object
      character(KIND=C_CHAR), intent(in)     :: gridName(*)
      integer(C_INT), intent(in), value      :: gridType
      integer(C_INT), intent(in), value      :: collectionType
    end function XdmfAddGrid

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Close the current opened XdmfGrid.
    !> @details
    !> Close the current open collection.  If within a nested collection, close
    !> the most deeply nested collection.
    !>
    !> Although this call is not mandatory because resources are freed when
    !> ~XdmfInterface is called, it is recommendable to call
    !> CloseGrid after WriteGrid and finishing to build the grid.
    !>
    !> @param[in]  object         Variable containing the XdmfInterface object ID
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfCloseGrid( object ) bind(c, name='XdmfCloseGrid')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfCloseGrid
      integer(C_INTPTR_T), intent(in), value :: object
    end function XdmfCloseGrid

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Set Topology associated with the current opened XdmfGrid.
    !> @details
    !> Sets the topology type to be assigned to the Grid.
    !>
    !> NumberType can be: \n
    !>        XDMF_INT_32 (INTEGER*4) \n
    !>        XDMF_INT_64 (INTEGER*8) \n
    !>
    !> TopologyType can be: \n
    !>        XDMF_NOTOPOLOGY     0x0 \n
    !>        XDMF_POLYVERTEX     0x1 \n
    !>        XDMF_POLYLINE       0x2 \n
    !>        XDMF_POLYGON        0x3 \n
    !>        XDMF_TRI            0x4 \n
    !>        XDMF_QUAD           0x5 \n
    !>        XDMF_TET            0x6 \n
    !>        XDMF_PYRAMID        0x7 \n
    !>        XDMF_WEDGE          0x8 \n
    !>        XDMF_HEX            0x9 \n
    !>        XDMF_EDGE_3         0x0022 \n
    !>        XDMF_TRI_6          0x0024 \n
    !>        XDMF_QUAD_8         0x0025 \n
    !>        XDMF_TET_10         0x0026 \n
    !>        XDMF_PYRAMID_13     0x0027 \n
    !>        XDMF_WEDGE_15       0x0028 \n
    !>        XDMF_WEDGE_18       0x0029 \n
    !>        XDMF_HEX_20         0x0030 \n
    !>        XDMF_HEX_24         0x0031 \n
    !>        XDMF_HEX_27         0x0032 \n
    !>        XDMF_MIXED          0x0070 \n
    !>        XDMF_2DSMESH        0x0100 \n
    !>        XDMF_2DRECTMESH     0x0101 \n
    !>        XDMF_2DCORECTMESH   0x0102 \n
    !>        XDMF_3DSMESH        0x1100 \n
    !>        XDMF_3DRECTMESH     0x1101 \n
    !>        XDMF_3DCORECTMESH   0x1102 \n
    !>
    !> Structured: \n
    !>  * 2DSMesh, 2DRectMesh, 2DCORECTMesh, 3DSMesh, 3DRectMesh and 3DCORECTMesh \n
    !>      - NumberOfElements: NumberOfElements (topology is implicit) \n
    !>
    !> Unstructured: \n
    !>  * POLY* topologies \n
    !>      - NodesPerElement:  NodesPerElement \n
    !>      - NumberOfElements: NumberOfElements \n
    !>      - Dimensions:       Dimensions (NodesPerElement * NumberOfElements) \n
    !>
    !>  * Mixed topologies \n
    !>      - NumberOfElements: NumberOfElements \n
    !>      - Dimensions:       Total number of nodes in mesh including elements descriptors \n
    !>
    !>  * Triangle, Quadrilateral, Tetrahedron, Pyramid, Wedge, Hexahedron, \n
    !>    Edge_3, Triangle_6, Quadrilateral_8, Tetrahedron_10, Pyramid_13, \n
    !>    Wedge_15, Wedge_18, Hexahedron_20, Hexahedron_24, Hexahedron_27 \n
    !>      - NumberOfElements: NumberOfElements (Dimensions/NodesPerElement) \n
    !>      - Dimensions:       Total number of nodes in mesh including elements descriptors \n
    !> 
    !> @param[in]    object       Variable containing the XdmfInterface object ID
    !> @param[inout] reference    Pointer to a reference ID object
    !> @param[in]    topologyType Type of Topology (see TopologyType description above)
    !> @param[in]    numberType   Type of numbers defining the Topology (see NumberType description above)
    !> @param[in]    noeRank      Rank of the numberOfElements variable
    !> @param[in]    numberOfElements Array of INT*8 containing the shape of the number of elements (elements)
    !> @param[in]    dimRank      Rank of the dimensions variable
    !> @param[in]    dimensions   Array of INT*8 containing the shape of the dimensions (nodes)
    !> @param[in]    conns        Array of numberType type containing the Topology (nodes)
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfSetGridTopology( object, reference, topologyType, numberType, noeRank, &
                                  numberOfElements, dimRank, dimensions, conns ) bind(c, name='XdmfSetGridTopology')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(conns)
      integer(C_INT)                         :: XdmfSetGridTopology
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      integer(C_INT), intent(in), value      :: topologyType
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: noeRank
      integer(C_LONG_LONG), intent(in)       :: numberOfElements
      integer(C_INT), intent(in), value      :: dimRank
      integer(C_LONG_LONG), intent(in)       :: dimensions
      integer(C_INTPTR_T), intent(in)        :: conns
    end function XdmfSetGridTopology

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Set Topology associated with the current opened XdmfGrid.
    !> @details
    !> Sets the topology type to be assigned to the Grid. From Shape definition version.
    !>
    !> @param[in]    object       Variable containing the XdmfInterface object ID
    !> @param[inout] reference    Pointer to a reference ID object
    !> @param[in]    topologyType Type of Topology (see TopologyType description above)
    !> @param[in]    numberType   Type of numbers defining the Topology (see NumberType description above)
    !> @param[in]    numberOfElements String containing the shape of the number of elements (elements)
    !> @param[in]    dimensions   String containing the shape of the dimensions (nodes)
    !> @param[in]    conns        Array of numberType type containing the Topology (nodes)
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfSetGridTopologyFromShape( object, reference, topologyType, numberType, &
                                           numberOfElements, dimensions, conns ) bind(c, name='XdmfSetGridTopologyFromShape')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(conns)
      integer(C_INT)                         :: XdmfSetGridTopologyFromShape
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      integer(C_INT), intent(in), value      :: topologyType
      integer(C_INT), intent(in), value      :: numberType
      character(KIND=C_CHAR), intent(in)     :: numberOfElements(*)
      character(KIND=C_CHAR), intent(in)     :: dimensions(*)
      integer(C_INTPTR_T), intent(in)        :: conns
    end function XdmfSetGridTopologyFromShape

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Set Geometry associated with the current opened XdmfGrid.
    !> @details
    !> Sets the geometry type to be assigned to the Grid.
    !>
    !> NumberType can be: \n
    !>        XDMF_FLOAT32 (REAL*4) \n
    !>        XDMF_FLOAT64 (REAL*8) \n
    !>
    !> TopologyType can be: \n
    !>        XDMF_GEOMETRY_NONE          0 \n
    !>        XDMF_GEOMETRY_XYZ           1 \n
    !>        XDMF_GEOMETRY_XY            2 \n
    !>        XDMF_GEOMETRY_X_Y_Z         3 \n
    !>        XDMF_GEOMETRY_X_Y           4 \n
    !>        XDMF_GEOMETRY_VXVYVZ        5 \n
    !>        XDMF_GEOMETRY_ORIGIN_DXDYDZ 6 \n
    !>        XDMF_GEOMETRY_VXVY          7 \n
    !>        XDMF_GEOMETRY_ORIGIN_DXDY   8 \n
    !>
    !> @param[in]    object       Variable containing the XdmfInterface object ID
    !> @param[inout] reference    Pointer to a reference ID object
    !> @param[in]    geometryType Type of Geometry (see TopologyType description above)
    !> @param[in]    numberType   Type of numbers defining the Topology (see NumberType description above)
    !> @param[in]    rank         Rank of the dimensions variable
    !> @param[in]    dimensions   Array of INT*8 containing the shape of the dimensions (nodes)
    !> @param[in]    points       Array of numberType type containing the Geometry (nodes)
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfSetGridGeometry( object, reference, geometryType, numberType, rank, &
                                  dimensions, points ) bind(c, name='XdmfSetGridGeometry')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(points)
      integer(C_INT)                         :: XdmfSetGridGeometry
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      integer(C_INT), intent(in), value      :: geometryType
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: rank
      integer(C_LONG_LONG), intent(in)       :: dimensions
      integer(C_INTPTR_T), intent(in)        :: points
    end function XdmfSetGridGeometry

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Set Geometry associated with the current opened XdmfGrid.
    !> @details
    !> Sets the geometry type to be assigned to the Grid. From Shape definition version.
    !>
    !> NumberType can be: \n
    !>        XDMF_FLOAT32 (REAL*4) \n
    !>        XDMF_FLOAT64 (REAL*8) \n
    !>
    !> TopologyType can be: \n
    !>        XDMF_GEOMETRY_NONE          0 \n
    !>        XDMF_GEOMETRY_XYZ           1 \n
    !>        XDMF_GEOMETRY_XY            2 \n
    !>        XDMF_GEOMETRY_X_Y_Z         3 \n
    !>        XDMF_GEOMETRY_X_Y           4 \n
    !>        XDMF_GEOMETRY_VXVYVZ        5 \n
    !>        XDMF_GEOMETRY_ORIGIN_DXDYDZ 6 \n
    !>        XDMF_GEOMETRY_VXVY          7 \n
    !>        XDMF_GEOMETRY_ORIGIN_DXDY   8 \n
    !>
    !> @param[in]    object       Variable containing the XdmfInterface object ID
    !> @param[inout] reference    Pointer to a reference ID object
    !> @param[in]    geometryType Type of Geometry (see GeometryType description above)
    !> @param[in]    numberType   Type of numbers defining the Geometry (see NumberType description above)
    !> @param[in]    dimensions   String containing the shape of the dimensions (nodes)
    !> @param[in]    points       Array of numberType type containing the Geometry (nodes)
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfSetGridGeometryFromShape( object, reference, geometryType, numberType, &
                                           dimensions, points ) bind(c, name='XdmfSetGridGeometryFromShape')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(points)
      integer(C_INT)                         :: XdmfSetGridGeometryFromShape
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      integer(C_INT), intent(in), value      :: geometryType
      integer(C_INT), intent(in), value      :: numberType
      character(KIND=C_CHAR), intent(in)     :: dimensions(*)
      integer(C_INTPTR_T), intent(in)        :: points
    end function XdmfSetGridGeometryFromShape

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Add attribute associated with the current opened XdmfGrid.
    !> @details
    !> Add an attribute to be written to the next opened grid. Multiple attributes can
    !> be added and written to a single grid.
    !>
    !> NumberTypes can be: \n
    !>         XDMF_INT8_TYPE      1 \n
    !>         XDMF_INT16_TYPE     6 \n
    !>         XDMF_INT32_TYPE     2 \n
    !>         XDMF_INT64_TYPE     3 \n
    !>         XDMF_FLOAT32_TYPE   4 \n
    !>         XDMF_FLOAT64_TYPE   5 \n
    !>         XDMF_UINT8_TYPE     7 \n
    !>         XDMF_UINT16_TYPE    8 \n
    !>         XDMF_UINT32_TYPE    9 \n
    !>
    !> AttributeCenter can be: \n
    !>         XDMF_ATTRIBUTE_CENTER_GRID  0 \n
    !>         XDMF_ATTRIBUTE_CENTER_CELL  1 \n
    !>         XDMF_ATTRIBUTE_CENTER_FACE  2 \n
    !>         XDMF_ATTRIBUTE_CENTER_EDGE  3 \n
    !>         XDMF_ATTRIBUTE_CENTER_NODE  4 \n
    !>
    !> AttributeType can be: \n
    !>         XDMF_ATTRIBUTE_TYPE_NONE       0x00 \n
    !>         XDMF_ATTRIBUTE_TYPE_SCALAR     0x01 \n
    !>         XDMF_ATTRIBUTE_TYPE_VECTOR     0x02 \n
    !>         XDMF_ATTRIBUTE_TYPE_TENSOR     0x03 \n
    !>         XDMF_ATTRIBUTE_TYPE_MATRIX     0x04 \n
    !>         XDMF_ATTRIBUTE_TYPE_TENSOR6    0x05 \n
    !>         XDMF_ATTRIBUTE_TYPE_GLOBALID   0x06 \n
    !>         XDMF_ATTRIBUTE_TYPE_UNIFORM    0x00 // By default Attributes are uniform \n
    !>         XDMF_ATTRIBUTE_TYPE_COLLECTION 0x10 // Attribute composed of several DataItems \n
    !>         XDMF_ATTRIBUTE_TYPE_MASK       0x0F // Evaluates type of Single DataItem \n
    !>
    !>    i.e.: XDMF_ATTRIBUTE_TYPE_SCALAR | XDMF_ATTRIBUTE_TYPE_COLLECTION
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[inout] reference       Pointer to a reference ID object
    !> @param[in]    attributeName   Name of the attribute
    !> @param[in]    numberType      Type of numbers defining the Attribute (see NumberType description above)
    !> @param[in]    attributeCenter Where the attribute is centered (see AttributeCenter description above)
    !> @param[in]    attributeType   Type of Attribute (see AttributeType description above)
    !> @param[in]    rank            Rank of the dimensions variable
    !> @param[in]    dimensions      Array of INT*8 containing the shape of the dimensions (nodes)
    !> @param[in]    data            Array of numberType type containing the Attribute
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfAddGridAttribute( object, reference, attributeName, numberType, attributeCenter, &
                                   attributeType, rank, dimensions, data ) bind(c, name='XdmfAddGridAttribute')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(data)
      integer(C_INT)                         :: XdmfAddGridAttribute
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      character(KIND=C_CHAR), intent(in)     :: attributeName(*)
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: attributeCenter
      integer(C_INT), intent(in), value      :: attributeType
      integer(C_INT), intent(in), value      :: rank
      integer(C_LONG_LONG), intent(in)       :: dimensions
      integer(C_INTPTR_T), intent(in)        :: data
    end function XdmfAddGridAttribute

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Add attribute associated with the current opened XdmfGrid.
    !> @details
    !> Add an attribute to be written to the next opened grid. Multiple attributes can
    !> be added and written to a single grid. From Shape description version.
    !>
    !> NumberTypes can be: \n
    !>         XDMF_INT8_TYPE      1 \n
    !>         XDMF_INT16_TYPE     6 \n
    !>         XDMF_INT32_TYPE     2 \n
    !>         XDMF_INT64_TYPE     3 \n
    !>         XDMF_FLOAT32_TYPE   4 \n
    !>         XDMF_FLOAT64_TYPE   5 \n
    !>         XDMF_UINT8_TYPE     7 \n
    !>         XDMF_UINT16_TYPE    8 \n
    !>         XDMF_UINT32_TYPE    9 \n
    !>
    !> AttributeCenter can be: \n
    !>         XDMF_ATTRIBUTE_CENTER_GRID  0 \n
    !>         XDMF_ATTRIBUTE_CENTER_CELL  1 \n
    !>         XDMF_ATTRIBUTE_CENTER_FACE  2 \n
    !>         XDMF_ATTRIBUTE_CENTER_EDGE  3 \n
    !>         XDMF_ATTRIBUTE_CENTER_NODE  4 \n
    !>
    !> AttributeType can be: \n
    !>         XDMF_ATTRIBUTE_TYPE_NONE       0x00 \n
    !>         XDMF_ATTRIBUTE_TYPE_SCALAR     0x01 \n
    !>         XDMF_ATTRIBUTE_TYPE_VECTOR     0x02 \n
    !>         XDMF_ATTRIBUTE_TYPE_TENSOR     0x03 \n
    !>         XDMF_ATTRIBUTE_TYPE_MATRIX     0x04 \n
    !>         XDMF_ATTRIBUTE_TYPE_TENSOR6    0x05 \n
    !>         XDMF_ATTRIBUTE_TYPE_GLOBALID   0x06 \n
    !>         XDMF_ATTRIBUTE_TYPE_UNIFORM    0x00 // By default Attributes are uniform \n
    !>         XDMF_ATTRIBUTE_TYPE_COLLECTION 0x10 // Attribute composed of several DataItems \n
    !>         XDMF_ATTRIBUTE_TYPE_MASK       0x0F // Evaluates type of Single DataItem \n
    !>
    !>    i.e.: XDMF_ATTRIBUTE_TYPE_SCALAR | XDMF_ATTRIBUTE_TYPE_COLLECTION
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[inout] reference       Pointer to a reference ID object
    !> @param[in]    attributeName   Name of the attribute
    !> @param[in]    numberType      Type of numbers defining the Attribute (see NumberType description above)
    !> @param[in]    attributeCenter Where the attribute is centered (see AttributeCenter description above)
    !> @param[in]    attributeType   Type of Attribute (see AttributeType description above)
    !> @param[in]    dimensions      String containing the shape of the dimensions (nodes)
    !> @param[in]    data            Array of numberType type containing the Attribute
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfAddGridAttributeFromShape( object, reference, attributeName, numberType, attributeCenter, &
                                            attributeType, shape, units, data ) bind(c, name='XdmfAddGridAttributeFromShape')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(data)
      integer(C_INT)                         :: XdmfAddGridAttributeFromShape
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      character(KIND=C_CHAR), intent(in)     :: attributeName(*)
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: attributeCenter
      integer(C_INT), intent(in), value      :: attributeType
      character(KIND=C_CHAR), intent(in)     :: shape(*)
      character(KIND=C_CHAR), intent(in)     :: units(*)
      integer(C_INTPTR_T), intent(in)        :: data
    end function XdmfAddGridAttributeFromShape

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Close the current opened XdmfAttribute.
    !> @details
    !> Close the current opened Attribute collection. If within a nested collection, close
    !> the most deeply nested collection.
    !>
    !> Although this call is not mandatory because resources are freed when
    !> ~XdmfInterface is called, it is recommendable to call
    !> CloseAttribute before closing Grid and WriteGrid.
    !>
    !> @param[in]  object         Variable containing the XdmfInterface object ID
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfCloseAttribute( object ) bind(c, name='XdmfCloseAttribute')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfCloseAttribute
      integer(C_INTPTR_T), intent(in), value :: object
    end function XdmfCloseAttribute

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Add a DataItem associated with the current opened XdmfAttribute.
    !> @details
    !> Add a DataItem to be written to the previous Attribute element. Multiple dataitems can
    !> be added and written into a single Attribute (Collection Attribute).
    !>
    !> NumberType can be: \n
    !>         XDMF_INT8_TYPE      1 \n
    !>         XDMF_INT16_TYPE     6 \n
    !>         XDMF_INT32_TYPE     2 \n
    !>         XDMF_INT64_TYPE     3 \n
    !>         XDMF_FLOAT32_TYPE   4 \n
    !>         XDMF_FLOAT64_TYPE   5 \n
    !>         XDMF_UINT8_TYPE     7 \n
    !>         XDMF_UINT16_TYPE    8 \n
    !>         XDMF_UINT32_TYPE    9 \n
    !>
    !> ItemType can be: \n
    !>         XDMF_ITEM_UNIFORM        0x00 \n
    !>         XDMF_ITEM_HYPERSLAB      0x01 \n
    !>         XDMF_ITEM_COORDINATES    0x02 \n
    !>         XDMF_ITEM_FUNCTION       0x03 \n
    !>         XDMF_ITEM_COLLECTION     0x14 \n
    !>         XDMF_ITEM_TREE           0x15 \n
    !>
    !>         XDMF_ITEM_MASK           0xF0    // Evaluates to a Single Array ? \n
    !>
    !> Format available - Only makes sense in XDMF_ITEM_UNIFORM: \n
    !>         XDMF_FORMAT_XML      0 \n
    !>         XDMF_FORMAT_HDF      1 \n
    !>         XDMF_FORMAT_MYSQL    2 \n
    !>         XDMF_FORMAT_BINARY   3 \n
    !>
    !> Reference is not returned back for a DataItem when externalHDF is enabled
    !> and it is a XDMF_ITEM_UNIFORM (Single Array) because there is not a real
    !> XdmfDataItem added to the XdmfInterface Tree actually.
    !> Future improvement? Add a real XdmfDataItem in this case?
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[inout] reference       Pointer to a reference ID object
    !> @param[in]    itemName        Name of the DataItem
    !> @param[in]    numberType      Type of numbers defining the DataItem (see NumberType description above)
    !> @param[in]    itemType        Type of DataItem (see ItemType description above)
    !> @param[in]    format          Specifies the output format of the DataItem (see Format description above)
    !> @param[in]    function        String that defines the function for a XDMF_ITEM_FUNCTION type
    !> @param[in]    rank            Rank of the dimensions variable
    !> @param[in]    dimensions      Array of INT*8 containing the shape of the dimensions
    !> @param[in]    data            Array of numberType type containing the DataItem element
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfAddDataItem( object, reference, itemName, numberType, itemType, &
                              format, function, rank, dimensions, data ) bind(c, name='XdmfAddDataItem')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(data)
      integer(C_INT)                         :: XdmfAddDataItem
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      character(KIND=C_CHAR), intent(in)     :: itemName(*)
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: itemType
      integer(C_INT), intent(in), value      :: format
      character(KIND=C_CHAR), intent(in)     :: function(*)
      integer(C_INT), intent(in), value      :: rank
      integer(C_LONG_LONG), intent(in)       :: dimensions
      integer(C_INTPTR_T), intent(in)        :: data
    end function XdmfAddDataItem

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Add a DataItem associated with the current opened XdmfAttribute.
    !> @details
    !> Add a DataItem to be written to the previous Attribute element. Multiple dataitems can
    !> be added and written into a single Attribute (Collection Attribute). From Shape definition version. 
    !>
    !> NumberType can be: \n
    !>         XDMF_INT8_TYPE      1 \n
    !>         XDMF_INT16_TYPE     6 \n
    !>         XDMF_INT32_TYPE     2 \n
    !>         XDMF_INT64_TYPE     3 \n
    !>         XDMF_FLOAT32_TYPE   4 \n
    !>         XDMF_FLOAT64_TYPE   5 \n
    !>         XDMF_UINT8_TYPE     7 \n
    !>         XDMF_UINT16_TYPE    8 \n
    !>         XDMF_UINT32_TYPE    9 \n
    !>
    !> ItemType can be: \n
    !>         XDMF_ITEM_UNIFORM        0x00 \n
    !>         XDMF_ITEM_HYPERSLAB      0x01 \n
    !>         XDMF_ITEM_COORDINATES    0x02 \n
    !>         XDMF_ITEM_FUNCTION       0x03 \n
    !>         XDMF_ITEM_COLLECTION     0x14 \n
    !>         XDMF_ITEM_TREE           0x15 \n
    !>
    !>         XDMF_ITEM_MASK           0xF0    // Evaluates to a Single Array ? \n
    !>
    !> Format available - Only makes sense in XDMF_ITEM_UNIFORM: \n
    !>         XDMF_FORMAT_XML      0 \n
    !>         XDMF_FORMAT_HDF      1 \n
    !>         XDMF_FORMAT_MYSQL    2 \n
    !>         XDMF_FORMAT_BINARY   3 \n
    !>
    !> Reference is not returned back for a DataItem when externalHDF is enabled
    !> and it is a XDMF_ITEM_UNIFORM (Single Array) because there is not a real
    !> XdmfDataItem added to the XdmfInterface Tree actually.
    !> Future improvement? Add a real XdmfDataItem in this case?
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[inout] reference       Pointer to a reference ID object
    !> @param[in]    itemName        Name of the DataItem
    !> @param[in]    numberType      Type of numbers defining the DataItem (see NumberType description above)
    !> @param[in]    itemType        Type of DataItem (see ItemType description above)
    !> @param[in]    format          Specifies the output format of the DataItem (see Format description above)
    !> @param[in]    function        String that defines the function for a XDMF_ITEM_FUNCTION type
    !> @param[in]    dimensions      String containing the shape of the dimensions
    !> @param[in]    data            Array of numberType type containing the DataItem element
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfAddDataItemFromShape( object, reference, itemName, numberType, itemType, &
                                       format, function, shape, data ) bind(c, name='XdmfAddDataItemFromShape')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(data)
      integer(C_INT)                         :: XdmfAddDataItemFromShape
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      character(KIND=C_CHAR), intent(in)     :: itemName(*)
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: itemType
      integer(C_INT), intent(in), value      :: format
      character(KIND=C_CHAR), intent(in)     :: function(*)
      character(KIND=C_CHAR), intent(in)     :: shape(*)
      integer(C_INTPTR_T), intent(in)        :: data
    end function XdmfAddDataItemFromShape

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Close the current opened XdmfDataItem.
    !> @details
    !> Close the current open DataItem collection. If within a nested collection, close
    !> the most deeply nested collection.
    !>
    !> Although this call is not mandatory because resources are freed when
    !> ~XdmfInterface is called, it is recommendable to call
    !> CloseDataItem before closing Grid and WriteGrid.
    !>
    !> @param[in]  object         Variable containing the XdmfInterface object ID
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfCloseDataItem( object ) bind(c, name='XdmfCloseDataItem')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfCloseDataItem
      integer(C_INTPTR_T), intent(in), value :: object
    end function XdmfCloseDataItem

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Add an Information element associated with the current opened Xdmf element.
    !> @details
    !> Add an Information element to be written to the last Grid/Attribute/DataItem element
    !> of the current collection.
    !> If we are not within a collection add to the top level domain
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[inout] reference       Pointer to a reference ID object
    !> @param[in]    informationName Name of the Information element
    !> @param[in]    value           String containing the information to be stored
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfAddInformation( object, reference, informationName, value ) bind(c, name='XdmfAddInformation')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfAddInformation
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      character(KIND=C_CHAR), intent(in)     :: informationName(*)
      character(KIND=C_CHAR), intent(in)     :: value(*)
    end function XdmfAddInformation

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Add a Set associated with the current opened XdmfGrid.
    !> @details
    !> Adds a XdmfSet element to the current opened Grid. Check XdmfSet element
    !> for further information about the structure and features.
    !>
    !> SetType (Where Ids are Assigned) can be: \n
    !>   XDMF_SET_TYPE_UNSET -1 \n
    !>   XDMF_SET_TYPE_NODE   1 \n
    !>   XDMF_SET_TYPE_CELL   2 \n
    !>      - Only one DataItem is written \n
    !> \n
    !>   XDMF_SET_TYPE_FACE   3 \n
    !>      - Two DataItems are written: \n
    !>          1st: Define CellIds \n
    !>          2nd: Define FaceIds in given CellIds \n
    !> \n
    !>   XDMF_SET_TYPE_EDGE   4 \n
    !>      - Three DataItems are written: \n
    !>          1st: Define CellIds \n
    !>          2nd: Define FaceIds in given CellIds \n
    !>          3rd: Define Edge in given FaceIds \n
    !>
    !>  When externalHDF is enabled several HDF5 URIs are separated through semicolons: \n
    !>    i.e.: test.h5:/CellIds;test.h5:/FaceIds;test.h5:/EdgeIds;test.h5:/Data
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[inout] reference       Pointer to a reference ID object
    !> @param[in]    setName         Name of the Set
    !> @param[in]    numberType      Type of numbers defining the Set
    !> @param[in]    setType         Type of Set (see SetType description above)
    !> @param[in]    attributeType   Type of Attribute for data (see AttributeType description above)
    !> @param[in]    rank            Rank of the dimensions variable
    !> @param[in]    dimensions      Array of INT*8 containing the shape of the dimensions
    !> @param[in]    cellIds         Array of numberType type containing the cell ids (XDMF_SET_TYPE_FACE, XDMF_SET_TYPE_EDGE)
    !> @param[in]    faceIds         Array of numberType type containing the face ids (XDMF_SET_TYPE_EDGE)
    !> @param[in]    ids             Array of numberType type containing the ids (cell ids, face ids, edge ids, node ids)
    !> @param[in]    data            Array of numberType type containing the DataItem element associated to the XdmfSet
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfAddGridSet( object, reference, setName, numberType, setType, attributeType, &
                             rank, dimensions, cellIds, faceIds, ids, data ) bind(c, name='XdmfAddGridSet')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(data)
      integer(C_INT)                         :: XdmfAddGridSet
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      character(KIND=C_CHAR), intent(in)     :: setName(*)
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: setType
      integer(C_INT), intent(in), value      :: attributeType
      integer(C_INT), intent(in), value      :: rank
      integer(C_LONG_LONG), intent(in)       :: dimensions
      integer(C_INTPTR_T), intent(in)        :: cellIds
      integer(C_INTPTR_T), intent(in)        :: faceIds
      integer(C_INTPTR_T), intent(in)        :: ids
      integer(C_INTPTR_T), intent(in)        :: data
    end function XdmfAddGridSet

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Add a Set associated with the current opened XdmfGrid.
    !> @details
    !> Adds a XdmfSet element to the current opened Grid. Check XdmfSet element
    !> for further information about the structure and features. From Shape definition version.
    !>
    !> SetType (Where Ids are Assigned) can be: \n
    !>   XDMF_SET_TYPE_UNSET -1 \n
    !>   XDMF_SET_TYPE_NODE   1 \n
    !>   XDMF_SET_TYPE_CELL   2 \n
    !>      - Only one DataItem is written \n
    !> \n
    !>   XDMF_SET_TYPE_FACE   3 \n
    !>      - Two DataItems are written: \n
    !>          1st: Define CellIds \n
    !>          2nd: Define FaceIds in given CellIds \n
    !> \n
    !>   XDMF_SET_TYPE_EDGE   4 \n
    !>      - Three DataItems are written: \n
    !>          1st: Define CellIds \n
    !>          2nd: Define FaceIds in given CellIds \n
    !>          3rd: Define Edge in given FaceIds \n
    !>
    !>  When externalHDF is enabled several HDF5 URIs are separated through semicolons: \n
    !>    i.e.: test.h5:/CellIds;test.h5:/FaceIds;test.h5:/EdgeIds;test.h5:/Data
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[inout] reference       Pointer to a reference ID object
    !> @param[in]    setName         Name of the Set
    !> @param[in]    numberType      Type of numbers defining the Set
    !> @param[in]    setType         Type of Set (see SetType description above)
    !> @param[in]    attributeType   Type of Attribute for data (see AttributeType description above)
    !> @param[in]    dimensions      String containing the shape of the dimensions
    !> @param[in]    cellIds         Array of numberType type containing the cell ids (XDMF_SET_TYPE_FACE, XDMF_SET_TYPE_EDGE)
    !> @param[in]    faceIds         Array of numberType type containing the face ids (XDMF_SET_TYPE_EDGE)
    !> @param[in]    ids             Array of numberType type containing the ids (cell ids, face ids, edge ids, node ids)
    !> @param[in]    data            Array of numberType type containing the DataItem element associated to the XdmfSet
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfAddGridSetFromShape( object, reference, setName, numberType, setType, attributeType, &
                                      shape, cellIds, faceIds, ids, data ) bind(c, name='XdmfAddGridSetFromShape')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(data)
      integer(C_INT)                         :: XdmfAddGridSetFromShape
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INTPTR_T), intent(inout)     :: reference
      character(KIND=C_CHAR), intent(in)     :: setName(*)
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: setType
      integer(C_INT), intent(in), value      :: attributeType
      character(KIND=C_CHAR), intent(in)     :: shape(*)
      integer(C_INTPTR_T), intent(in)        :: cellIds
      integer(C_INTPTR_T), intent(in)        :: faceIds
      integer(C_INTPTR_T), intent(in)        :: ids
      integer(C_INTPTR_T), intent(in)        :: data
    end function XdmfAddGridSetFromShape

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Read a Xdmf file.
    !> @details
    !> Reads a Xdmf file into the current XdmfDOM.  Must call XdmfReadGrid() to read
    !> in associated geometry, topology, and attributes.
    !> To be TESTED!
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[in]    filePath        Name of the Xdmf file (.xmf) to be read
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfReadFile( object, filePath ) bind(c, name='XdmfReadFile')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
      character(KIND=C_CHAR), intent(in)     :: filePath(*)
    end subroutine XdmfReadFile

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Read a XdmfGrid.
    !> @details
    !> Reads a Grid in the current XdmfDOM into XdmfGeometry, XdmfTopology,
    !> and XdmfAttribute elements. Grid selected by gridName.
    !> An XdmfReadGrid() followed by a XdmfWriteGrid() will make a copy of the grid.
    !> To be TESTED!
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[in]    gridName        Name of the XdmfGrid element to read
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfReadGrid( object, gridName ) bind(c, name='XdmfReadGrid')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
      character(KIND=C_CHAR), intent(in)     :: gridName(*)
    end subroutine XdmfReadGrid

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Read a XdmfGrid.
    !> @details
    !> Reads a Grid in the current XdmfDOM into XdmfGeometry, XdmfTopology,
    !> and XdmfAttribute elements. Grid is selected by index order in XdmfDOM.
    !> An XdmfReadGrid() followed by a XdmfWriteGrid() will make a copy of the grid.
    !> To be TESTED!
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[in]    gridIndex       Index of the XdmfGrid element to read
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfReadGridAtIndex( object, gridIndex ) bind(c, name='XdmfReadGridAtIndex')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INT), intent(in), value      :: gridIndex
    end subroutine XdmfReadGridAtIndex

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Get the Number of Grids in the current XdmfDOM.
    !> @detail
    !> Returns the number of grids in the current open file. This ignores collections.
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @return Number of Grids in the current XdmfDOM
    !---------------------------------------------------------------------------  
    function XdmfGetNumberOfGrids( object ) bind(c, name='XdmfGetNumberOfGrids')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfGetNumberOfGrids
      integer(C_INTPTR_T), intent(in), value :: object
    end function XdmfGetNumberOfGrids

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Get the Number of Elements in the current open XdmfGrid.
    !> @detail
    !> Returns the number of elements in the current open Grid (the current active
    !> XdmfTopology Element). This is either from a current read-in file or from
    !> a created but unwritten grid. If no topology element is present, return -1.
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[out]   dimensions      Array of INT*8 containing the shape of the Number Of Elements
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfGetNumberOfElements( object, dimensions ) bind(c, name='XdmfGetNumberOfElements')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfGetNumberOfElements
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_LONG_LONG), intent(out)      :: dimensions
    end function XdmfGetNumberOfElements

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Get the Number of Points in the current open XdmfGrid.
    !> @detail
    !> Returns the number of points in the current open grid (the current active
    !> XdmfGeometry Element). This is either from a current read-in file or from
    !> a created but unwritten grid.  If no geometry element is present, return -1.
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @return Number of Points in the current active XdmfGeometry Element
    !---------------------------------------------------------------------------  
    function XdmfGetNumberOfPoints( object ) bind(c, name='XdmfGetNumberOfPoints')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfGetNumberOfPoints
      integer(C_INTPTR_T), intent(in), value :: object
    end function XdmfGetNumberOfPoints

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Read the point values from the current XdmfGeometry.
    !> @detail
    !> Reads the point values from the current geometry into the passed array
    !> pointer.  If the geometry has not been created no values are read.
    !> To be TESTED!
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[in]    numberType      Type of numbers defining the Point Values
    !> @param[in]    startIndex      Starting point of the array to fill
    !> @param[out]   arrayToFill     Array of numberType type to fill with point values from XdmfGeometry
    !> @param[in]    numberOfValues  Number of elements to read
    !> @param[in]    arrayStride     Stride for arrayToFill buffer
    !> @param[in]    valueStride     Stride for point values in XdmfGeometry
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfReadPointValues( object, numberType, startIndex, arrayToFill, &
                                    numberOfValues, arrayStride, valuesStride ) bind(c, name='XdmfReadPointValues')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(arrayToFill)
      integer(C_INTPTR_T), intent(in), value :: object
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: startIndex
      integer(C_INTPTR_T), intent(out)       :: arrayToFill
      integer(C_INT), intent(in), value      :: numberOfValues
      integer(C_INT), intent(in), value      :: arrayStride
      integer(C_INT), intent(in), value      :: valuesStride
    end subroutine XdmfReadPointValues

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Get the number of values in the specified XdmfAttribute.
    !> @detail
    !> Returns the number of values in the specified attribute. Iterates over all
    !> current open attributes to find the specified attribute name and returns the
    !> number of values it contains.  If no attribute is found, return -1.
    !> To be TESTED!
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[in]    attributeName   String containing the name of the XdmfAttribute
    !> @return Number of values for attributeName XdmfAttribute
    !---------------------------------------------------------------------------  
    function XdmfGetNumberOfAttributeValues( object, attributeName ) bind(c, name='XdmfGetNumberOfAttributeValues')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfGetNumberOfAttributeValues
      integer(C_INTPTR_T), intent(in), value :: object
      character(KIND=C_CHAR), intent(in)     :: attributeName(*)
    end function XdmfGetNumberOfAttributeValues

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Read the values from the specified XdmfAttribute.
    !> @detail
    !> Reads the values from the specified attribute into the passed array pointer.
    !> If the attribute cannot be found, no values are read.
    !> To be TESTED!
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[in]    attributeName   String containing the name of the XdmfAttribute
    !> @param[in]    numberType      Type of numbers defining the Attribute Values
    !> @param[in]    startIndex      Starting point of the array to fill
    !> @param[out]   arrayToFill     Array of numberType type to fill with values from XdmfAttribute
    !> @param[in]    numberOfValues  Number of elements to read
    !> @param[in]    arrayStride     Stride for arrayToFill buffer
    !> @param[in]    valueStride     Stride for point values in XdmfAttribute
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfReadAttributeValues( object, attributeName, numberType, startIndex, arrayToFill, &
                                        numberOfValues, arrayStride, valuesStride ) bind(c, name='XdmfReadAttributeValues')
      use iso_c_binding
      implicit none
      !NO_ARG_CHECK(arrayToFill)
      integer(C_INTPTR_T), intent(in), value :: object
      character(KIND=C_CHAR), intent(in)     :: attributeName(*)
      integer(C_INT), intent(in), value      :: numberType
      integer(C_INT), intent(in), value      :: startIndex
      integer(C_INTPTR_T), intent(out)       :: arrayToFill
      integer(C_INT), intent(in), value      :: numberOfValues
      integer(C_INT), intent(in), value      :: arrayStride
      integer(C_INT), intent(in), value      :: valuesStride
    end subroutine XdmfReadAttributeValues

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Read the values from the specified XdmfInformation element.
    !> @detail
    !> Return the values from the specified information element.  If the
    !> information element cannot be found, no values are passed.  Information
    !> elements at the top level domain are searched first, followed by the
    !> currently loaded grid.
    !> To be TESTED!
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[in]    informationName String containing the name of the XdmfInformation
    !> @return String of the XdmfInformation element
    !---------------------------------------------------------------------------  
    function XdmfReadInformationValue( object, informationName ) bind(c, name='XdmfReadInformationValue')
      use iso_c_binding
      implicit none
      character(KIND=C_CHAR)                 :: XdmfReadInformationValue
      integer(C_INTPTR_T), intent(in), value :: object
      character(KIND=C_CHAR), intent(in)     :: informationName(*)
    end function XdmfReadInformationValue

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Builds XdmfDOM structure with the current definition of the Grid.
    !> @detail
    !> Add a grid to the XdmfDOM.  Assign the current topology, geometry, and grid
    !> attributes to grid.  If within a collection, add grid to the collection,
    !> otherwise add to the top level domain.  Assign time value if value is
    !> nonnegative.
    !> This function must be call after the definition of each Grid in order
    !> to generate the Xml datastructure for a following output.
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @return status error (XDMF_SUCCESS or XDMF_FAIL)
    !---------------------------------------------------------------------------  
    function XdmfWriteGrid( object ) bind(c, name='XdmfWriteGrid')
      use iso_c_binding
      implicit none
      integer(C_INT)                         :: XdmfWriteGrid
      integer(C_INTPTR_T), intent(in), value :: object
    end function XdmfWriteGrid

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Writes constructed Xdmf file (.xmf) to disk.
    !> @detail
    !> Dumps built Xml Xdmf memory structure to disk with filename created
    !> upon initialization. This function should be called after XdmfWriteGrid.
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfWriteToFile( object ) bind(c, name='XdmfWriteToFile')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
    end subroutine XdmfWriteToFile

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Serialize XdmfDOM structure to stdout.
    !> @detail
    !> Prints current XdmfDOM to console
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfSerialize( object ) bind(c, name='XdmfSerialize')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
    end subroutine XdmfSerialize

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Serialize XdmfDOM structure to stdout.
    !> @detail
    !> Copy current XdmfDOM object to memory pointed by charPointer
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @param[inout] charPointer     Pointer to the address where copying XdmfDOM object ID
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfGetDOM( object, charPointer ) bind(c, name='XdmfGetDOM')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
      character(C_CHAR), intent(out)         :: charPointer(*)
    end subroutine XdmfGetDOM

    !---------------------------------------------------------------------------  
    ! DESCRIPTION: 
    !> @brief Close XdmfInterface interface and clean memory.
    !> @detail
    !> Closes and destroys the XdmfInterface cleaning and freeing the memory.
    !> This function must be called at the end of the Xdmf execution.
    !>
    !> @param[in]    object          Variable containing the XdmfInterface object ID
    !> @return nothing
    !---------------------------------------------------------------------------  
    subroutine XdmfClose( object ) bind(c, name='XdmfClose')
      use iso_c_binding
      implicit none
      integer(C_INTPTR_T), intent(in), value :: object
    end subroutine XdmfClose
  end interface

  contains

  !---------------------------------------------------------------------------  
  ! DESCRIPTION: 
  !> @brief Converts a Fortran 90 string to C format string.
  !> @detail
  !> Helper function to convert form Fortran-format strings to C-format.
  !>
  !> @param[in]    F_STRING   String in Fortran90 format
  !> @return C_STRING, String in C format
  !---------------------------------------------------------------------------  
  pure function C_STR(F_STRING) result(C_STRING)
    use, intrinsic :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
    implicit none
    character(LEN=*), intent(in) :: F_STRING
    character(LEN=1,KIND=C_CHAR) :: C_STRING(LEN_TRIM(F_STRING)+1)
    INTEGER                      :: N, I

    N = LEN_TRIM(F_STRING)
    DO I = 1, N
      C_STRING(I) = F_STRING(I:I)
    END DO
    C_STRING(N + 1) = C_NULL_CHAR

    !C_STRING = F_STRING // C_NULL_CHAR
  end function C_STR

end module Xdmf
