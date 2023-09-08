!------------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @name    Coupling functions
!> @file    mod_couplings.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for coupli
!> @details ToolBox for coupli 
!>          To create a coupling:
!>
!>          1. Initialize the structure:
!>          - call COU_INITIALIZE_COUPLING_STRUCTURE(COUPLING)
!>          2. Compute the following:
!>          - COUPLING % GEOME % NUMBER_WET_POINTS ....................... Number of wet points
!>          - COUPLING % GEOME % NPOIN_WET ............................... Number of wet nodes
!>          - COUPLING % GEOME % COORD_WET(:,1:NUMBER_WET_POINTS) ........ Coordinates of wet points
!>          - COUPLING % GEOME % LPOIN_WET(1:NPOIN_WET) .................. List of wet nodes
!>          - COUPLING % ITYPE ........................................... Vector projection
!>          - COUPLING % KIND ............................................ BETWEEN_SUBDOMAINS/BETWEEN_ZONES
!>          - COUPLING % COLOR_TARGET .................................... Target color
!>          - COUPLING % COLOR_SOURCE .................................... Source color
!>          3. Initialize the coupling:
!>          - call COU_INIT_INTERPOLATE_POINTS_VALUES(coupling_type(icoup) % geome % coord_wet,color_target,color_source,COUPLING)
!>          4. Generate transmission matrices:
!>          - call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling,what_to_do)
!>          - call COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES(coupling,what_to_do)
!>          - call COU_PARALLELIZE_TRANSMISSION_MATRICES(coupling)
!> 
!>          Trick, if the NUMBER_WET_POINTS wet points do not have anything to do with the
!>          mesh, do the following:
!>
!>          - KIND = BETWEEN_ZONES
!>          - NPOIN_WET = NUMBER_WET_POINTS
!>          - LPOIN_WET(1:NPOIN_WET) = 1:NPOIN_WET
!>
!>          Implicit coupling
!>          -----------------
!>
!>          * Mass matrix
!>
!>          Let x = [1,1,1...]^t
!>          The lumped mass matrix Ml can be obtained from the mass matrix M as
!>          Ml = M.x
!>          Thus it can be computed just like the classical SpMV y = Ax
!>          which is valid as far as xd=Td.xn. This is our case a 1d=Td.1n due
!>          to the constant conservation of the Dirichlet transmission matrix.
!>           
!> @{
!------------------------------------------------------------------------

module mod_couplings
  use def_kintyp,         only : ip,rp,lg,r1p,r2p,i1p,i2p,spmat,comm_data_par  
  use def_master,         only : ISEQUEN
  use def_master,         only : intost
  use def_master,         only : ittim
  use def_master,         only : current_code
  use def_master,         only : current_zone
  use def_master,         only : INOTMASTER 
  use def_master,         only : IMASTER,kfl_paral,lninv_loc
  use def_master,         only : zeror,lzone,namda
  use def_master,         only : I_AM_IN_SUBD
  use def_domain,         only : lnods
  use def_domain,         only : lesub,nelem,lnnod
  use def_domain,         only : mnodb,mnode,ndime
  use def_domain,         only : nbono,coord,ltopo,ltype
  use def_domain,         only : lbono,npoin,lnoch
  use def_domain,         only : nboun,lnodb,lnnob
  use def_domain,         only : lelch,meshe
  use def_domain,         only : nnode
  use def_kermod,         only : ielse,relse,ndivi
  use def_elmtyp,         only : NOFRI
  use def_elmtyp,         only : ELHOL
  use mod_elmgeo,         only : elmgeo_natural_coordinates
  use mod_elmgeo,         only : elmgeo_natural_coordinates_on_boundaries
  use mod_elsest,         only : elsest_host_element
  use mod_parall,         only : PAR_THIS_NODE_IS_MINE
  use mod_parall,         only : PAR_COMM_COLOR_PERM
  use mod_parall,         only : par_part_in_color
  use mod_parall,         only : par_code_zone_subd_to_color
  use mod_parall,         only : color_target
  use mod_parall,         only : color_source
  use mod_parall,         only : par_bin_comin
  use mod_parall,         only : par_bin_comax
  use mod_parall,         only : par_bin_part
  use mod_parall,         only : par_bin_boxes
  use mod_parall,         only : par_bin_size
  use mod_parall,         only : PAR_COMM_CURRENT
  use mod_parall,         only : PAR_COMM_COLOR
  use mod_parall,         only : PAR_COMM_COLOR_ARRAY
  use mod_parall,         only : I_AM_IN_COLOR
  use mod_parall,         only : PAR_MY_CODE_RANK
  use mod_parall,         only : par_part_comin
  use mod_parall,         only : par_part_comax
  use mod_parall,         only : PAR_GLOBAL_TO_LOCAL_NODE
  use mod_parall,         only : PAR_INITIALIZE_COMMUNICATION_ARRAY
  use mod_parall,         only : PAR_WORLD_SIZE
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_alloca_min
  use mod_memory,         only : memory_size
  use mod_memory,         only : memory_resize
  use mod_interpolation,  only : COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_interpolation,  only : COU_GENERATE_PROJECTION_MATRIX
  use mod_kdtree,         only : typ_kdtree
  use mod_kdtree,         only : kdtree_nearest_boundary
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_BARRIER
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_START_NON_BLOCKING_COMM
  use mod_communications, only : PAR_END_NON_BLOCKING_COMM
  use mod_communications, only : PAR_SET_NON_BLOCKING_COMM_NUMBER
  use mod_communications, only : PAR_GATHER
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_ALLGATHERV
  use mod_communications, only : PAR_SEND_RECEIVE_TO_ALL
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only : PAR_ALLTOALL
  use def_coupli,         only : kdtree_typ
  use def_coupli,         only : mcoup
  use def_coupli,         only : coupling_type
  use def_coupli,         only : UNKNOWN
  use def_coupli,         only : RESIDUAL
  use def_coupli,         only : DIRICHLET_IMPLICIT
  use def_coupli,         only : DIRICHLET_EXPLICIT 
  use def_master,         only : AT_BEGINNING  
  use def_coupli,         only : RELAXATION_SCHEME
  use def_coupli,         only : UNKNOWN
  use def_coupli,         only : ELEMENT_INTERPOLATION
  use def_coupli,         only : BOUNDARY_INTERPOLATION
  use def_coupli,         only : NEAREST_BOUNDARY_NODE
  use def_coupli,         only : NEAREST_ELEMENT_NODE
  use def_coupli,         only : BOUNDARY_VECTOR_PROJECTION
  use def_coupli,         only : GLOBAL_NUMBERING
  use def_coupli,         only : SAME_COORDINATE  
  use def_coupli,         only : TRANSPOSE_MIRROR
  use def_coupli,         only : ON_WHOLE_MESH
  use def_coupli,         only : typ_color_coupling
  use def_coupli,         only : memor_cou
  use def_coupli,         only : RELAXATION_SCHEME
  use def_coupli,         only : AITKEN_SCHEME
  use def_coupli,         only : BROYDEN_SCHEME
  use def_coupli,         only : IQNLS_SCHEME
  use def_coupli,         only : STRESS_PROJECTION
  use def_coupli,         only : PROJECTION
  use def_coupli,         only : BETWEEN_SUBDOMAINS
  use def_coupli,         only : BETWEEN_ZONES
  use def_coupli,         only : ON_CHIMERA_MESH
  use def_coupli,         only : FIXED_UNKNOWN
  use def_coupli,         only : INTERFACE_MASS  
  use def_coupli,         only : GLOBAL_MASS    
  use def_coupli,         only : ncoup_implicit 
  use def_coupli,         only : coupling_driver_couplings
  use def_coupli,         only : coupling_driver_iteration
  use def_coupli,         only : coupling_driver_max_iteration
  use def_coupli,         only : coupling_driver_number_couplings
  use def_coupli,         only : coupling_driver_tolerance
  use def_coupli,         only : nboun_cou
  use def_coupli,         only : lnodb_cou
  use def_coupli,         only : ltypb_cou
  use def_coupli,         only : lboch_cou
  use def_coupli,         only : lnnob_cou
  use mod_matrix,         only : matrix_COO_spgemm
  use mod_matrix,         only : matrix_COO_aggregate
  use mod_matrix,         only : nullify_spmat
  use mod_iofile,         only : iofile_open_unit
  use mod_iofile,         only : iofile_close_unit
  use mod_iofile,         only : iofile_available_unit
  use mod_iofile,         only : iofile_flush_unit
  use mod_messages,       only : messages_live
  use mod_maths,          only : maths_matrix_vector_multiplication
  use mod_maths,          only : maths_outer_product
  use mod_std
  implicit none 
  private
 
!!$  interface COU_INTERPOLATE_NODAL_VALUES
!!$     module procedure COU_INTERPOLATE_NODAL_VALUES_11_real,&
!!$          &           COU_INTERPOLATE_NODAL_VALUES_22_real,&
!!$          &           COU_INTERPOLATE_NODAL_VALUES_33_real,&
!!$          &           COU_INTERPOLATE_NODAL_VALUES_12_real,&
!!$          &           COU_INTERPOLATE_NODAL_VALUES_21_real,&
!!$          &           COU_INTERPOLATE_NODAL_VALUES_1_pointer,&
!!$          &           COU_INTERPOLATE_NODAL_VALUES_2_pointer
!!$  end interface COU_INTERPOLATE_NODAL_VALUES

  integer(ip), parameter :: my_huge = huge(1_ip)
  
  interface COU_PUT_VALUE_ON_TARGET
     module procedure COU_PUT_VALUE_ON_TARGET_IP_1,&  
          &           COU_PUT_VALUE_ON_TARGET_IP_2,&
          &           COU_PUT_VALUE_ON_TARGET_IP_12
  end interface COU_PUT_VALUE_ON_TARGET
  
  interface COU_INIT_INTERPOLATE_POINTS_VALUES
     module procedure COU_INIT_INTERPOLATE_POINTS_VALUES_OLD,&
          &           COU_INIT_INTERPOLATE_POINTS_VALUES_NEW
  end interface COU_INIT_INTERPOLATE_POINTS_VALUES

  public :: COU_INTERPOLATE_NODAL_VALUES
  public :: COU_RESIDUAL_FORCE
  public :: I_AM_IN_COUPLING
  public :: I_AM_INVOLVED_IN_A_COUPLING_TYPE
  public :: COU_INIT_INTERPOLATE_POINTS_VALUES                ! Initialize color coupling
  public :: COU_PRESCRIBE_DIRICHLET_IN_MATRIX
  public :: COU_LIST_SOURCE_NODES
  public :: COU_CHECK_CONVERGENCE
  public :: THERE_EXISTS_A_ZONE_COUPLING
  public :: I_HAVE_A_FRINGE_ELEMENT
  public :: MIRROR_COUPLING
  public :: I_HAVE_A_CHIMERA_COUPLING
  public :: COU_PUT_VALUE_ON_TARGET
  public :: COU_SET_FIXITY_ON_TARGET
  public :: COU_TEMPORAL_PREDICTOR
  public :: couplings_initialize_solver
  public :: couplings_impose_dirichlet
  public :: couplings_check_dirichlet                         ! Check if an array satisfies the Dirichlet coupling condition
  public :: couplings_initialization
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-26
  !> @brief   Couplingsinitialization
  !> @details Initializaiton of coupling variables
  !> 
  !-----------------------------------------------------------------------

  subroutine couplings_initialization()

    use def_coupli
    !
    ! Memory
    !
    memor_cou              = 0
    !
    ! Read variables: initialization must be here as sometimes
    ! couplign is used without data file
    !
    mcoup                  = 0
    toler_absolute_cou     = -1.00_rp                ! Element average size for partition bounding box
    toler_relative_cou     =  0.01_rp                ! Relative tolerance for partition bounding box
    kfl_lost_wet_point_cou = STOP_ALYA_WITH_WARNINGS ! Stop Alya with warnings
    number_of_holes        = 0                       ! No hole
    !
    ! Driver
    !  
    coupling_driver_couplings     = 0
    coupling_driver_max_iteration = 1      
    coupling_driver_iteration     = 0   
    coupling_driver_tolerance     = 1.0e-12_rp    
    !
    ! Others
    !
    nullify(coupling_type)
    !
    ! Implicit couplings
    !
    ncoup_implicit   = 0
    ncoup_implicit_d = 0   
    ncoup_implicit_n = 0   
    nullify(lcoup_implicit_d) 
    nullify(lcoup_implicit_n) 
    nullify(mask_cou)      

  end subroutine couplings_initialization
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    02/10/2014
  !> @brief   Find mirror coupling
  !> @details Obtain the mirror coupling MIRROR_COUPLING of ICOUP
  !>          It returns zero if no mirror has been found
  !>
  !----------------------------------------------------------------------

  function I_HAVE_A_CHIMERA_COUPLING()
    integer(ip) :: icoup 
    logical(lg) :: I_HAVE_A_CHIMERA_COUPLING

    I_HAVE_A_CHIMERA_COUPLING = .false.
    do icoup = 1,mcoup
       if(    I_AM_IN_COUPLING(icoup) .and. &
            & coupling_type(icoup) % where_type == ON_CHIMERA_MESH ) then
          I_HAVE_A_CHIMERA_COUPLING = .true.
          return
       end if
    end do

  end function I_HAVE_A_CHIMERA_COUPLING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    02/10/2014
  !> @brief   Find mirror coupling
  !> @details Obtain the mirror coupling MIRROR_COUPLING of ICOUP
  !>          It returns zero if no mirror has been found
  !
  !----------------------------------------------------------------------

  function MIRROR_COUPLING(icoup)
    integer(ip), intent(in) :: icoup        !< Coupling
    integer(ip)             :: MIRROR_COUPLING
    logical(lg)             :: notfound

    MIRROR_COUPLING = 0
    notfound = .true.
    do while( notfound .and. MIRROR_COUPLING < mcoup )
       MIRROR_COUPLING = MIRROR_COUPLING + 1
       if(  &
            & coupling_type(icoup) % color_target           == coupling_type(MIRROR_COUPLING) % color_source .and. &
            & coupling_type(MIRROR_COUPLING) % color_target == coupling_type(icoup) % color_source ) then
          coupling_type(icoup) % mirror_coupling = MIRROR_COUPLING
          notfound = .false.
       end if
    end do

  end function MIRROR_COUPLING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   If I am in a coupling
  !> @details Check if I am involved in cupling icoup
  !
  !----------------------------------------------------------------------

  function I_AM_IN_COUPLING(icoup)
    integer(ip), intent(in) :: icoup        !< Coupling
    integer(ip)             :: icolo_source
    integer(ip)             :: icolo_target
    logical(lg)             :: I_AM_IN_COUPLING

    icolo_source = coupling_type(icoup) % color_source
    icolo_target = coupling_type(icoup) % color_target
    I_AM_IN_COUPLING = I_AM_IN_COLOR(icolo_source) .or. I_AM_IN_COLOR(icolo_target)

  end function I_AM_IN_COUPLING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   If I am in a coupling of a certain type
  !> @details If I am in a coupling of a certain type: RESIDUAL, UNKNOWN
  !>          DIRICHLET_EXPLICIT or DIRICHLET_IMPLICIT
  !
  !----------------------------------------------------------------------

  function I_AM_INVOLVED_IN_A_COUPLING_TYPE(ikind,iwhat,itype,micou_type)
     integer(ip), intent(in)           :: ikind    !< Coupling kind (between zones, between subdomain)
     integer(ip), intent(in)           :: iwhat    !< Coupling what (residual, unknown, dirichlet)
     integer(ip), intent(in), optional :: itype    !< Coupling type (projection, interpolation)
     integer(ip), intent(in), optional :: micou_type !< Mirror coupling type (projection, interpolation)
     integer(ip)                       :: icoup
     logical(lg)                       :: I_AM_INVOLVED_IN_A_COUPLING_TYPE
     integer(ip)                       :: micou

     I_AM_INVOLVED_IN_A_COUPLING_TYPE = .false.

     if( ikind == BETWEEN_SUBDOMAINS .and. ncoup_implicit == 0 ) then
        continue
     else
        do icoup = 1,mcoup
           if(    I_AM_IN_COUPLING(icoup) .and. &
              & coupling_type(icoup) % kind == ikind .and. &
              & coupling_type(icoup) % what == iwhat ) then
              I_AM_INVOLVED_IN_A_COUPLING_TYPE = .true.
              if(present(itype) .and. (.not.(coupling_type(icoup) % itype == itype))) &
                 & I_AM_INVOLVED_IN_A_COUPLING_TYPE = .false.
              if(present(micou_type)) then
                 micou = coupling_type(icoup) % mirror_coupling
                 if(.not.(coupling_type(micou) % itype == micou_type)) &
                    & I_AM_INVOLVED_IN_A_COUPLING_TYPE = .false.
              end if
           end if
        end do
     end if

  end function I_AM_INVOLVED_IN_A_COUPLING_TYPE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    23/09/2014
  !> @brief   Check if I have a fringe element
  !> @details A fringe element is a hole element with at least one if
  !>          fringe node
  !
  !----------------------------------------------------------------------

  function I_HAVE_A_FRINGE_ELEMENT()
    integer(ip) :: ipoin,ielem,inode
    logical(lg) :: I_HAVE_A_FRINGE_ELEMENT

    I_HAVE_A_FRINGE_ELEMENT = .false.
    do ielem = 1,nelem
       if( lelch(ielem) == ELHOL ) then
          do inode = 1,lnnod(ielem)
             ipoin = lnods(inode,ielem)
             if( lnoch(ipoin) == NOFRI ) then
                I_HAVE_A_FRINGE_ELEMENT = .true.
                return
             end if
          end do
       end if
    end do

  end function I_HAVE_A_FRINGE_ELEMENT

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    23/09/2014
  !> @brief   Check if I have a hole element
  !> @details A hole element has only hole nodes
  !
  !----------------------------------------------------------------------

  function I_HAVE_A_HOLE_ELEMENT()
    integer(ip) :: ipoin,ielem,inode
    logical(lg) :: I_HAVE_A_HOLE_ELEMENT

    I_HAVE_A_HOLE_ELEMENT = .false.
    do ielem = 1,nelem
       if( lelch(ielem) == ELHOL ) then
          I_HAVE_A_HOLE_ELEMENT = .true.
          return
       end if
    end do

  end function I_HAVE_A_HOLE_ELEMENT

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    30/09/2014
  !> @brief   Check convergence
  !> @details Check convergence of a coupling
  !
  !----------------------------------------------------------------------

  subroutine COU_CHECK_CONVERGENCE(iblok,kfl_gozon)
    integer(ip), intent(in)  :: iblok
    integer(ip), intent(out) :: kfl_gozon
    integer(ip)              :: icoup,kcoup
    integer(ip)              :: idime,idofn,ndofn,ipoin,kpoin
    integer(ip)              :: ntime_con,itime

    kfl_gozon = 0

    if( coupling_driver_iteration(iblok) >= coupling_driver_max_iteration(iblok) ) then
       !
       ! Not converged, max numb. of iteration overpassed
       !
       do kcoup = 1,coupling_driver_number_couplings(iblok)
          icoup = coupling_driver_couplings(kcoup,iblok)
          !
          ! Update values for subcycling coupling (frequency of exchanges different of one)
          !
          if( coupling_type(icoup) % frequ_send > 1_ip .or. coupling_type(icoup) % frequ_recv > 1_ip )then
             ! 
             ! Only source code save exchanging values
             !
             if( current_code == coupling_type(icoup) % code_source )then

                do kpoin = 1_ip, coupling_type(icoup) % geome % npoin_source
                   ipoin = coupling_type(icoup) % geome % lpoin_source(kpoin)
                   do idime = 1_ip, ndime
                      coupling_type(icoup) % values_frequ(idime,kpoin,1_ip) = coupling_type(icoup) % values_frequ(idime,kpoin,2_ip)
                   end do
                end do
             end if
          end if
          ! 
          ! Save last exchanged values for temporal predictor. Only the target code will make predictions
          !
          if( coupling_type(icoup) % kind == BETWEEN_ZONES .and. current_code == coupling_type(icoup) % code_target .and. coupling_type(icoup) % temporal_predictor == 1_ip )then
             ndofn = size(coupling_type(icoup) % values_converged,1_ip)
             ntime_con = size(coupling_type(icoup) % values_converged,3)
             do itime = ntime_con,2,-1

                do ipoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                   do idofn = 1,ndofn
                      coupling_type(icoup) % values_converged(idofn,ipoin,itime) = coupling_type(icoup) % values_converged(idofn,ipoin,itime-1)
                   end do
                end do

             end do
             do ipoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                do idofn = 1,ndofn
                   coupling_type(icoup) % values_converged(idofn,ipoin,1_ip) = coupling_type(icoup) % values(idofn,ipoin,1_ip)
                end do
             end do

          end if

       end do
       return
    else

       do kcoup = 1,coupling_driver_number_couplings(iblok)
          icoup = coupling_driver_couplings(kcoup,iblok)
          !
          ! Only check if the current time step is a multiple of the frequency defined
          !
          if( mod( ittim,coupling_type(icoup) % frequ_send ) == 0_ip .and. current_code == coupling_type(icoup) % code_source )then
             if( coupling_type(icoup) % resid(1) > coupling_driver_tolerance(iblok) )then 
                kfl_gozon = 1
             end if
          end if
          if( mod( ittim,coupling_type(icoup) % frequ_recv ) == 0_ip .and. current_code == coupling_type(icoup) % code_target )then
             if( coupling_type(icoup) % resid(1) > coupling_driver_tolerance(iblok) )then 
                kfl_gozon = 1
             end if
          end if
       end do
       !
       ! Update values exchanged when convergence is reached and frequency of exchanges is larger than one
       !
       do kcoup = 1,coupling_driver_number_couplings(iblok)
          icoup = coupling_driver_couplings(kcoup,iblok)

          if( coupling_type(icoup) % frequ_send > 1_ip .or. coupling_type(icoup) % frequ_recv > 1_ip )then
             ! do kcoup = 1,coupling_driver_number_couplings(iblok)
             !    icoup = coupling_driver_couplings(kcoup,iblok)
             if( mod( ittim,coupling_type(icoup) % frequ_send) == 0_ip .and. kfl_gozon==0 .and. current_code == coupling_type(icoup) % code_source )then
                do kpoin = 1_ip, coupling_type(icoup) % geome % npoin_source
                   ipoin = coupling_type(icoup) % geome % lpoin_source(kpoin)
                   do idime = 1_ip, ndime
                      coupling_type(icoup) % values_frequ(idime,kpoin,1_ip) = coupling_type(icoup) % values_frequ(idime,kpoin,2_ip)
                   end do
                end do
             end if
             ! end do
          end if





          ! 
          ! Save converged values for temporal predictor. Only the target code will make predictions
          !


          if( coupling_type(icoup) % kind == BETWEEN_ZONES .and. current_code == coupling_type(icoup) % code_target .and. kfl_gozon==0 .and. &
               coupling_type(icoup) % temporal_predictor == 1_ip )then

             ndofn = size(coupling_type(icoup) % values_converged,1_ip)
             ntime_con = size(coupling_type(icoup) % values_converged,3)


             ! do kcoup = 1,coupling_driver_number_couplings(iblok)
             !    icoup = coupling_driver_couplings(kcoup,iblok)

             do itime = ntime_con,2,-1
                do ipoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                   do idofn = 1,ndime
                      coupling_type(icoup) % values_converged(idofn,ipoin,itime) = coupling_type(icoup) % values_converged(idofn,ipoin,itime-1)
                   end do
                end do
             end do

             do ipoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                do idofn = 1,ndime
                   coupling_type(icoup) % values_converged(idofn,ipoin,1_ip) = coupling_type(icoup) % values(idofn,ipoin,1_ip)
                end do
             end do

             ! end do
          end if  ! coupling kind
       end do     ! loop over couplings
    end if


  end subroutine COU_CHECK_CONVERGENCE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   Actualize unknown according to scheme
  !> @details Actualize the unknown and save previous values according 
  !>          to the scheme (relaxation, Aitken, etc.)
  !
  !----------------------------------------------------------------------
  subroutine COU_UPDATE_POINTS_VALUES(xxnew,coupling,xresi,mask)
    use mod_iqnls

    implicit none 

    real(rp),      pointer,   intent(inout)        :: xxnew(:,:)
    type(typ_color_coupling), intent(inout)        :: coupling
    real(rp),                 intent(out)          :: xresi(2)
    integer(ip),   pointer,   intent(in), optional :: mask(:,:)
    integer(ip)                                    :: ndofn,ipoin,npoin_wet,ntime_wet
    integer(ip)                                    :: idofn,itime,kpoin,jpoin
    real(rp)                                       :: relax,rela1,rip1,rip2,xfact
    real(rp)                                       :: numer,denom,normv,vj,di,ti

    integer(ip)                                    :: actual_iter, i_iter, computation_iter, ndofs, idofs, mindex

    real(rp)                                       :: alpha(coupling % ranku_iqnls)
    real(rp)                                       :: reaux
    real(rp)                                       :: vecaux(coupling % wet % npoin_wet * size(xxnew,1) )
    
    coupling % itera = coupling % itera + 1
    ntime_wet = 0

    if( INOTMASTER ) then
       !
       ! Degrees of freedom
       !
       ndofn     = size(xxnew,1)
       npoin_wet = coupling % wet % npoin_wet
       ndofs      = npoin_wet * ndofn
       !
       ! Allocate memory if required and initializate some values
       !
       if( .not. associated( coupling % values_converged ) .and. npoin_wet > 0 .and. coupling  % temporal_predictor == 1_ip )then
          !
          ! Three time steps values will be stored to make temporal predictions the allocation is done here because it is not 
          ! know a priori the ndofn of the exchanged quantity
          !
          call memory_alloca(memor_cou,'MOD_COUPLINGS ',' values_converged ',coupling % values_converged,ndofn,npoin_wet,3_ip)
          coupling % values_converged = 0.0_rp
       end if

       if( coupling % scheme == RELAXATION_SCHEME )then     

          if( .not. associated(coupling % values) )then
             call memory_alloca(memor_cou,'MOD_COUPLINGS','values',coupling % values,ndofn,npoin_wet,2_ip)
          endif

          if( associated(coupling % values) ) ntime_wet = size(coupling % values,3,kind=ip)

       else if( coupling % scheme == AITKEN_SCHEME )then     

          if( .not. associated(coupling % values) )then
             call memory_alloca(memor_cou,'MOD_COUPLINGS','values',coupling % values,ndofn,npoin_wet,3_ip)
             call memory_alloca(memor_cou,'MOD_COUPLINGS','values_predicted',coupling % values_predicted,ndofn,npoin_wet)
          endif

          if( associated(coupling % values) ) ntime_wet = size(coupling % values,3,kind=ip)

       else if( coupling % scheme == BROYDEN_SCHEME )then     

          if( .not. associated(coupling % values) )then
             !             call memory_alloca(memor_cou,'MOD_COUPLINGS','values',coupling % values,ndofn,npoin_wet,3_ip)
             !             call memory_alloca(memor_cou,'MOD_COUPLINGS','values_predicted',coupling % values_predicted,ndofn,npoin_wet)
             !             call memory_alloca(memor_cou,'MOD_COUPLINGS','jacobian_inverse',coupling % jacobian_inverse,ndofn,npoin_wet,npoin_wet)
             !             call memory_alloca(memor_cou,'MOD_COUPLINGS','dincr_predicted',coupling % dincr_predicted,npoin_wet)
          endif

          if( associated(coupling % values) ) ntime_wet = size(coupling % values,3,kind=ip)

       else if( coupling % scheme == IQNLS_SCHEME )then     

          if( .not. associated(coupling % relaxed_iqnls) )then
             if(ndofs>0_ip) then
                call memory_alloca(memor_cou,'MOD_COUPLINGS','values',coupling % values,ndofn,npoin_wet,1_ip)

                call memory_alloca(memor_cou,'MOD_COUPLINGS','relaxed', coupling % relaxed_iqnls,           ndofn * npoin_wet, coupling % ranku_iqnls)
                call memory_alloca(memor_cou,'MOD_COUPLINGS','unrelaxed_iqnls', coupling % unrelaxed_iqnls, ndofn * npoin_wet, coupling % ranku_iqnls)
                call memory_alloca(memor_cou,'MOD_COUPLINGS','residues_iqnls', coupling % residues_iqnls,   ndofn * npoin_wet, coupling % ranku_iqnls)
                call memory_alloca(memor_cou,'MOD_COUPLINGS','valincr_iqnls',  coupling % valincr_iqnls,    ndofn * npoin_wet, coupling % ranku_iqnls)
                call memory_alloca(memor_cou,'MOD_COUPLINGS','residincr_iqnls', coupling % residincr_iqnls, ndofn * npoin_wet, coupling % ranku_iqnls)
             else
                call memory_alloca_min( coupling % values )
                call memory_alloca_min( coupling % values_converged )
                call memory_alloca_min( coupling % values_predicted )
                call memory_alloca_min( coupling % relaxed_iqnls )
                call memory_alloca_min( coupling % unrelaxed_iqnls )
                call memory_alloca_min( coupling % residues_iqnls )
                call memory_alloca_min( coupling % valincr_iqnls )
                call memory_alloca_min( coupling % residincr_iqnls )
             endif
          endif

          if( associated(coupling % relaxed_iqnls) ) ntime_wet = size(coupling % relaxed_iqnls,2)

       end if

    else
       ndofn     = 0_ip
       npoin_wet = 0_ip
       ntime_wet = 0_ip
       ndofs     = 0_ip

!!! ACAAAAAAAAAAAAAAAAA !TODO
!!! ACAAAAAAAAAAAAAAAAA !TODO
!!! ACAAAAAAAAAAAAAAAAA !TODO
       !!       if( coupling % scheme == IQNLS_SCHEME )then     
       call memory_alloca_min( coupling % values )
       call memory_alloca_min( coupling % values_converged )
       call memory_alloca_min( coupling % values_predicted )
       call memory_alloca_min( coupling % relaxed_iqnls )
       call memory_alloca_min( coupling % unrelaxed_iqnls )
       call memory_alloca_min( coupling % residues_iqnls )
       call memory_alloca_min( coupling % valincr_iqnls )
       call memory_alloca_min( coupling % residincr_iqnls )
!!! ACAAAAAAAAAAAAAAAAA !TODO
!!! ACAAAAAAAAAAAAAAAAA !TODO
!!! ACAAAAAAAAAAAAAAAAA !TODO
       !!       elseif( (coupling % scheme == AITKEN_SCHEME) .or. (coupling % scheme == RELAXATION_SCHEME) )then
       !!         call memory_alloca_min( coupling % values )
       !!         call memory_alloca_min( coupling % values_converged )
       !!         call memory_alloca_min( coupling % values_predicted )
       !!       endif



    end if


    !
    ! Save old relaxed values
    !
    ! relaxed(1)=x^k+1         <- to be calculated
    ! relaxed(2)=x^k        
    ! relaxed(3)=x^k-1
    ! relaxed(4)=x^k-2
    !        .
    !        .
    !        .
    !
    !
    ! And save old residues, if IQNLS
    ! residues_iqnls(1)=r^k    <- to be written
    ! residues_iqnls(2)=r^k-1
    ! residues_iqnls(3)=r^k-2
    !        .
    !        .
    !        .

!!! ACAAAAAAAAAAAAAAAAA !TODO
!!! ACAAAAAAAAAAAAAAAAA !TODO
!!! ACAAAAAAAAAAAAAAAAA !TODO
!!!  if( coupling % scheme .ne. BROYDEN_SCHEME) then

    if( coupling % scheme .eq. IQNLS_SCHEME) then
       !
       ! if IQNLS save also residues and unrelaxed values
       !
       do itime = ntime_wet,2,-1
          do idofs=1, ndofs
             coupling % relaxed_iqnls(idofs,itime)    = coupling % relaxed_iqnls(idofs,itime-1)
             coupling % unrelaxed_iqnls(idofs, itime) = coupling % unrelaxed_iqnls(idofs, itime-1)
             coupling % residues_iqnls(idofs, itime)  = coupling % residues_iqnls(idofs, itime-1)
          enddo
       end do

    elseif ( (coupling % scheme == RELAXATION_SCHEME) .or. (coupling % scheme == AITKEN_SCHEME)) then
       !
       ! Rest of the schemes but broyden
       !
       do itime = ntime_wet,2,-1
          do ipoin = 1,npoin_wet
             do idofn = 1,ndofn
                coupling % values(idofn,ipoin,itime) = coupling % values(idofn,ipoin,itime-1)
             end do
          end do
       end do

    endif
!!! ACAAAAAAAAAAAAAAAAA !TODO
!!! ACAAAAAAAAAAAAAAAAA !TODO
!!! ACAAAAAAAAAAAAAAAAA !TODO
!!!!  endif

    !
    ! Aitken: copy predicted value
    ! 
    if( coupling % scheme == RELAXATION_SCHEME ) then
       relax = coupling % relax
    else if( coupling % scheme == AITKEN_SCHEME ) then
       !
       ! First two iterations are performed with constant relaxation
       !
       if( coupling_driver_iteration(1_ip) < 3_ip ) then
          ! if( coupling % itera < 3_ip ) then

          coupling % aitken = coupling % relax
          relax             = coupling % aitken
          ! print*, "DEBUG: AITKEN inicial ", relax

       else
          !
          ! Scalar product for aitken relaxation factor
          !
          numer = 0.0_rp
          denom = 0.0_rp
          do ipoin = 1,npoin_wet
             do idofn = 1,ndofn
                !
                ! rip1 = d_{i}   - d_{i+1}'
                !
                rip1  = coupling % values_predicted(idofn,ipoin) - coupling % values(idofn,ipoin,3_ip) 
                ! rip1 = ( coupling % values(idofn,ipoin,2_ip) -  coupling % values(idofn,ipoin,3_ip) ) * xxnew(idofn,ipoin)
                !
                ! rip2 = d_{i+1} - d_{i+2}'
                !
                rip2  =  xxnew(idofn,ipoin) - coupling % values(idofn,ipoin,2_ip)
                ! rip2 = xxnew(idofn,ipoin) - coupling % values_predicted(idofn,ipoin)
                ! rip3 = xxnew(idofn,ipoin) - coupling % values(idofn,ipoin,2_ip)
                numer = numer + rip1 * (rip2-rip1)
                denom = denom + (rip2-rip1)*(rip2-rip1)
                ! numer = numer + rip1
                ! denom = denom + rip2 * rip3
             end do
          end do

          call PAR_SUM(numer,'IN CURRENT COLOR')
          call PAR_SUM(denom,'IN CURRENT COLOR')
          relax =-coupling % aitken * numer / (denom+zeror)
          ! if( relax < 0_rp .or. relax > 0.6 ) relax = 0.1_rp
          ! relax = - numer / (denom + zeror)
          ! relax = min(relax, 1._rp)
          ! relax = max(relax,-1._rp)
          coupling % aitken = relax
          ! if( relax < 0.01_rp ) relax = 0.02_rp
          ! if( relax > 1.00_rp ) relax = 0.03_rp
          ! aux = dabs(relax)
          ! print*, "DEBUG: AITKEN ", relax
       end if

    else if( coupling % scheme == BROYDEN_SCHEME ) then

       call COU_BROYDEN_BAD(xxnew,coupling)

    else if( coupling % scheme == IQNLS_SCHEME ) then

       actual_iter = coupling_driver_iteration(1_ip)

       !!
       !! Initialize values in first iteration
       !!
       if( actual_iter .eq. 1_ip ) then
          coupling % relaxed_iqnls=0.0_rp
          coupling % unrelaxed_iqnls=0.0_rp
          coupling % residues_iqnls=0.0_rp
          coupling % valincr_iqnls=0.0_rp
          coupling % residincr_iqnls=0.0_rp
          alpha = 0.0_rp
       endif

       !!
       !! Save new values and compute residues
       !!

       mindex=0_ip
       do ipoin=1, npoin_wet
          do idofn=1, ndofn
             mindex = ipoin * ndofn + idofn - ndofn

             ! In the first iteration of the current time step save the last
             ! converged value in the actual relaxed vector to compute
             ! the residue more precisely
             !
             if (actual_iter .eq. 1_ip .and. associated(coupling % values_converged) ) coupling % relaxed_iqnls(mindex,2) = coupling % values_converged(idofn,ipoin,1_ip)

             ! Save new result, for notation consistency
             coupling % unrelaxed_iqnls(mindex,1) = xxnew(idofn, ipoin) * coupling % scaling_iqnls
             ! Compute new residue
             coupling % residues_iqnls(mindex,1) = coupling % unrelaxed_iqnls(mindex,1) - coupling % relaxed_iqnls (mindex, 2)
          enddo
       enddo

       !
       ! First two iterations are performed with constant relaxation
       !
       if( actual_iter < 2_ip ) then

          relax = coupling % relax

          do idofs=1,ndofs
             coupling % relaxed_iqnls (idofs, 1) = relax * coupling % unrelaxed_iqnls(idofs, 1) + ( 1.0_rp -relax ) * coupling % relaxed_iqnls(idofs,2)
          enddo

       else
          !-----------------------------------------------------
          ! From the iteration 3, the actual IQN-LS 
          ! will be computed
          !
          !-----------------------------------------------------

          ! If we have more iterations than columns, use the
          ! maximum number of columns instead of the iterations
          !
          if (actual_iter .le. coupling % ranku_iqnls) then
             computation_iter = actual_iter
          else
             computation_iter = coupling % ranku_iqnls
          endif


          ! Build increment vector for this iteration
          !
          ! Build vector V_i= r_i - r_k
          ! and vector   W_i = x*_i - x*_k
          !
          do i_iter=1,computation_iter-1
             do idofs=1, ndofs
                coupling % residincr_iqnls(idofs, i_iter) = coupling % residues_iqnls(idofs,computation_iter-i_iter)  - coupling % residues_iqnls(idofs,computation_iter)
                coupling % valincr_iqnls(idofs, i_iter)   = coupling % unrelaxed_iqnls(idofs,computation_iter-i_iter) -  coupling % unrelaxed_iqnls(idofs,computation_iter)
             enddo
          enddo

          alpha = 0.0_rp

          call compute_alpha( coupling % residincr_iqnls(:,:) , &
               coupling % residues_iqnls(:,1)  , &
               ndofs, computation_iter-1, alpha )

          ! Matrix vector product W*alpha
          ! used for the new prediction
          !
          vecaux=0.0_rp
          call maths_matrix_vector_multiplication( coupling % valincr_iqnls(:,:), alpha(:), vecaux(:), computation_iter - 1 )

          ! New prediction
          ! x^k+1 = x^k + W*alpha -r^k
          !
          do idofs=1,ndofs
             coupling % relaxed_iqnls (idofs, 1) = coupling % relaxed_iqnls (idofs, 2) + coupling % residues_iqnls(idofs, 1) + vecaux(idofs)
          enddo

          !! Here finishes the differentiation between iter<3 and iter>3
       endif

       !! We must save again the result in xxnew
       !!
       do ipoin=1, npoin_wet
          do idofn=1, ndofn
             mindex = ipoin * ndofn + idofn - ndofn
             xxnew(idofn,ipoin) = coupling % relaxed_iqnls (mindex, 1) / coupling % scaling_iqnls
             coupling % values(idofn,ipoin,1) = xxnew(idofn,ipoin) 
          enddo
       enddo

       !! Here  finishes the SCHEME_IQNLS
    end if  ! RELAXATION SCHEME

    !
    ! Save unrelaxed results for next aitken calculation (must be performed in all iterations)
    !
    if(( coupling % scheme == AITKEN_SCHEME)) then
       do ipoin = 1_ip, npoin_wet
          do idofn = 1_ip, ndofn
             coupling % values_predicted(idofn,ipoin) = xxnew(idofn,ipoin)
          end do
       end do
    end if
    !
    ! Relaxation of the solution
    !    
    if( coupling % scheme == AITKEN_SCHEME .and. coupling_driver_iteration(1_ip) > 2_ip ) then
       rela1    = 1.0_rp - relax
       do kpoin = 1,npoin_wet
          do idofn = 1,ndofn
             xxnew(idofn,kpoin)               = rela1 * coupling % values(idofn,kpoin,2_ip) + relax * xxnew(idofn,kpoin)
             coupling % values(idofn,kpoin,1) = xxnew(idofn,kpoin)
          end do
       end do
    else if( coupling % scheme == BROYDEN_SCHEME ) then


    else if ( (coupling % scheme == RELAXATION_SCHEME) .or. &
         (coupling % scheme == AITKEN_SCHEME .and. coupling_driver_iteration(1_ip) <= 2_ip ) ) then
       rela1    = 1.0_rp - relax
       do kpoin = 1,npoin_wet
          ipoin = coupling % wet % lpoin_wet(kpoin)
          do idofn = 1,ndofn
             xxnew(idofn,kpoin)               = relax * xxnew(idofn,kpoin) + rela1 * coupling % values(idofn,kpoin,2)
             coupling % values(idofn,kpoin,1) = xxnew(idofn,kpoin)
          end do
       end do
    end if
    ! 
    ! Compute residual
    !
    xresi = 0.0_rp
    if( coupling % scheme == IQNLS_SCHEME ) then
       do kpoin = 1,npoin_wet
          ipoin = coupling % wet % lpoin_wet(kpoin)
          if( PAR_THIS_NODE_IS_MINE(ipoin) ) then
             do idofn = 1,ndofn
                mindex = kpoin * ndofn + idofn - ndofn
                xresi(2) = xresi(2) +  coupling % relaxed_iqnls(mindex,1)**2
                xresi(1) = xresi(1) + (coupling % relaxed_iqnls(mindex,1)-coupling % relaxed_iqnls(mindex,2))**2
             end do
          end if
       end do
    else
       xfact = 1.0_rp
       do kpoin = 1,npoin_wet
          ipoin = coupling % wet % lpoin_wet(kpoin)
          if( PAR_THIS_NODE_IS_MINE(ipoin) ) then
             do idofn = 1,ndofn
                if( present(mask) ) then
                   xfact = 1.0_rp
                   if(       coupling % what == RESIDUAL ) then
                      if( mask(idofn,ipoin) > 0 ) xfact = 0.0_rp
                   else if( coupling %  what == UNKNOWN ) then
                      if( mask(idofn,ipoin) > 0 .and. mask(idofn,ipoin) /= FIXED_UNKNOWN )  xfact = 0.0_rp
                  end if
                end if
                xresi(2) = xresi(2) + xfact *  coupling % values(idofn,kpoin,1)**2
                xresi(1) = xresi(1) + xfact * (coupling % values(idofn,kpoin,1)-coupling % values(idofn,kpoin,2))**2                   
             end do
          end if
       end do
    end if

  end subroutine COU_UPDATE_POINTS_VALUES

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Interpolate and modify nodal arrays
  !> @details Do the following
  !>                               Interpolate
  !>          XTARGET(NDOFN,NPOIN)     <=       XSOURCE(NDOFN,NPOIN)
  !>        
  !>          1. Allocate XINTERP(NDOFN,NPOIN_WET) 
  !>
  !>          2. Interpolate XINTERP(NDOFN,NPOIN_WET) FROM XSOURCE:
  !>
  !>                     COU_GET_INTERPOLATE_POINTS_VALUES 
  !>             XINTERP              <=                   XSOURCE
  !>             
  !>
  !>          3. Scatter solution:
  !>                     LPOIN_WET
  !>             XTARGET    <=    XINTERP
  !>
  !----------------------------------------------------------------------

  subroutine COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,xtarget,xsource,mask)

    use mod_parall, only : PAR_MY_WORLD_RANK
    
    integer(ip), intent(in)                    :: icoup
    integer(ip), intent(in)                    :: ndofn
    real(rp),    intent(inout)                 :: xtarget(ndofn,*) 
    real(rp),    intent(in),          optional :: xsource(ndofn,*) 
    integer(ip), intent(in), pointer, optional :: mask(:,:)
    real(rp),                pointer           :: xinterp(:,:)
    real(rp)                                   :: xresi(2)
    integer(ip)                                :: ipoin,kpoin,idofn,icolo,npoin_wet
    integer(ip)                                :: jcoup
    logical(lg)                                :: subdo_coupling
    real(rp)                                   :: weight,xnew,xold
    real(rp),    pointer                       :: xsource_loc(:,:)
    logical(lg)                                :: if_mask
    !
    ! Initialization
    !
    subdo_coupling = .false.
    if( coupling_type(icoup) % zone_target == 0 ) then
       if( I_AM_IN_SUBD( coupling_type(icoup) % subdomain_target) ) subdo_coupling = .true.
       if( I_AM_IN_SUBD( coupling_type(icoup) % subdomain_source) ) subdo_coupling = .true.
       icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)     
    else
       icolo = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
    end if

    color_target = coupling_type(icoup) % color_target
    color_source = coupling_type(icoup) % color_source
    npoin_wet    = coupling_type(icoup) % wet % npoin_wet
    xresi        = 0.0_rp
    if( IMASTER ) npoin_wet = 0
    nullify(xinterp)
    !
    ! Options
    !
    if_mask = .false.
    if( present(mask) ) then
       if( associated(mask) ) if_mask = .true.
    end if
    
    if( icolo == color_target .or. icolo == color_source .or. subdo_coupling ) then
       !
       ! Allocate values
       !
       if( current_code == coupling_type(icoup) % code_target .and. INOTMASTER ) then  
          call memory_alloca(memor_cou,'XINTERP','cou_interpolate_nodal_values',xinterp,ndofn,max(1_ip,npoin_wet))
       else
          call memory_alloca(memor_cou,'XINTERP','cou_interpolate_nodal_values',xinterp,1_ip,1_ip)
       end if
       !
       ! Interpolate values from XSOURCE
       !
       if( present(xsource) ) then
          call COU_GET_INTERPOLATE_POINTS_VALUES(ndofn,xsource,xinterp,coupling_type(icoup))
       else
          nullify(xsource_loc)
          call memory_alloca(memor_cou,'XSOURCE_LOC','cou_interpolate_nodal_values',xsource_loc,ndofn,max(1_ip,npoin))
          if( INOTMASTER ) xsource_loc(1:ndofn,1:npoin) = xtarget(1:ndofn,1:npoin)
          call COU_GET_INTERPOLATE_POINTS_VALUES(ndofn,xsource_loc,xinterp,coupling_type(icoup))
          call memory_deallo(memor_cou,'XSOURCE_LOC','cou_interpolate_nodal_values',xsource_loc)
       end if
       !
       ! Permute value to target XTARGET
       !
       if( current_code == coupling_type(icoup) % code_target ) then
          !
          ! Update solution according to scheme (relaxation, Aitken, etc.)
          !
          if( .not. subdo_coupling ) call COU_UPDATE_POINTS_VALUES(xinterp,coupling_type(icoup),xresi,mask)

          if(    coupling_type(icoup) %  what == UNKNOWN            .or. &
               & coupling_type(icoup) %  what == DIRICHLET_IMPLICIT .or. &
               & coupling_type(icoup) %  what == DIRICHLET_EXPLICIT ) then

             if( if_mask ) then
                !
                ! Unknown with mask
                !
                do kpoin = 1,npoin_wet
                   ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                   do idofn = 1,ndofn
                      if( mask(idofn,ipoin) <= 0 .or. mask(idofn,ipoin) == FIXED_UNKNOWN ) then
                         xtarget(idofn,ipoin) = xinterp(idofn,kpoin)
                      end if
                   end do
                end do
             else
                !
                ! Unknown without mask
                !
                do kpoin = 1,npoin_wet
                   ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                   do idofn = 1,ndofn
                      xtarget(idofn,ipoin) = xinterp(idofn,kpoin)
                   end do
                end do
             end if
             
          else if( coupling_type(icoup) % what == RESIDUAL ) then

             if( if_mask ) then
                !
                ! Force with mask
                !
                do kpoin = 1,npoin_wet
                   ipoin  = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                   weight = coupling_type(icoup) % wet % weight_wet(kpoin)
                   do idofn = 1,ndofn 
                      if( mask(idofn,ipoin) <= 0 ) then
                         xtarget(idofn,ipoin) = xtarget(idofn,ipoin) + weight * xinterp(idofn,kpoin)
                      end if
                   end do
                end do
             else
                !
                ! Force without mask
                !
                do kpoin = 1,npoin_wet
                   ipoin  = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                   weight = coupling_type(icoup) % wet % weight_wet(kpoin)
                   do idofn = 1,ndofn
                     xtarget(idofn,ipoin) = xtarget(idofn,ipoin) + weight * xinterp(idofn,kpoin)
                   end do
                end do
             end if
          else
             call runend('WRONG TAG')
          end if
          !
          ! Conservation
          !
          if( if_mask ) then
             call COU_CONSERVATION(coupling_type(icoup),ndofn,xtarget,mask)
          end if
       end if
       call memory_deallo(memor_cou,'XINTERP','cou_interpolate_nodal_values',xinterp)
       !
       ! It should be in color icolo and jcolo!
       !
       if( .not. subdo_coupling ) then
          call PAR_SUM(2_ip,xresi,'IN CURRENT COUPLING')
          coupling_type(icoup) % resid(2) = sqrt(xresi(2))
          coupling_type(icoup) % resid(1) = sqrt(xresi(1)) / ( sqrt(xresi(2)) + zeror )
       end if 

    end if

  end subroutine COU_INTERPOLATE_NODAL_VALUES

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Conservation
  !> @details Apply a conservation algorithm
  !>
  !----------------------------------------------------------------------

  subroutine COU_CONSERVATION(coupling,ndofn,xtarget,mask)
    use mod_projec, only : projec_mass_conservation
    type(typ_color_coupling), intent(in)    :: coupling
    integer(ip),              intent(in)    :: ndofn
    integer(ip),              intent(in)    :: mask(ndofn,*)
    real(rp),                 intent(inout) :: xtarget(ndofn,*) 
    integer(ip)                             :: kboun,iboun,inodb,ipoin
    integer(ip)                             :: idofn
    logical(lg)                             :: mark_node
    logical(lg),              pointer       :: gboun(:)
    logical(lg),              pointer       :: gpoin(:)

    if( coupling % conservation /= 0 ) then

       if( INOTMASTER ) then
          allocate( gboun(nboun) )
          allocate( gpoin(npoin) )
          gboun = .false.
          gpoin = .false.
          do kboun = 1,coupling % wet % nboun_wet 
             iboun = coupling % wet % lboun_wet(kboun)
             gboun(iboun) = .true.
             do inodb = 1,lnnob(iboun)
                ipoin = lnodb(inodb,iboun) 
                mark_node = .true.
                loop_idofn: do idofn = 1,ndofn
                   if( mask(idofn,ipoin) /= FIXED_UNKNOWN ) then
                      mark_node = .false.
                      exit loop_idofn
                   end if
                end do loop_idofn
                if( mark_node ) gpoin(ipoin) = .true.
             end do
          end do
       else
          allocate( gboun(1) )
          allocate( gpoin(1) )
       end if
       if( coupling % conservation == INTERFACE_MASS ) then
          call projec_mass_conservation(xtarget,gboun,gpoin,'LOCAL MASS')
       else if( coupling % conservation == GLOBAL_MASS ) then
          call projec_mass_conservation(xtarget,gboun,gpoin,'GLOBAL MASS')    
       else
          call runend('COU_CONSERVATION: CONSERVATION NOT CODED')
       end if
       deallocate( gboun )
       deallocate( gpoin )

    end if

  end subroutine COU_CONSERVATION

  subroutine COU_RESIDUAL_FORCE(ndofn,ia,ja,amatr,rhsid,unkno,force_rhs)
    use def_kintyp, only        :  ip,rp
    use def_master, only        :  INOTMASTER
    use def_domain, only        :  npoin
    implicit none
    integer(ip),    intent(in)  :: ndofn
    integer(ip),    intent(in)  :: ia(*)
    integer(ip),    intent(in)  :: ja(*)
    real(rp),       intent(in)  :: amatr(ndofn,ndofn,*)
    real(rp),       intent(in)  :: rhsid(ndofn,*)
    real(rp),       intent(in)  :: unkno(ndofn,*)
    real(rp),       intent(out) :: force_rhs(ndofn,*)
    integer(ip)                 :: ipoin,iz,jpoin,idofn,jdofn
    !
    ! Marcar unicamente los nodos que tocan los elementos
    !
    if( INOTMASTER ) then
       do ipoin = 1,npoin
          force_rhs(1:ndofn,ipoin) = rhsid(1:ndofn,ipoin)
          do iz = ia(ipoin),ia(ipoin+1)-1
             jpoin = ja(iz)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   force_rhs(idofn,ipoin) = &
                        force_rhs(idofn,ipoin) - amatr(jdofn,idofn,iz)*unkno(jdofn,jpoin)
                end do
             end do
          end do
       end do
       call rhsmod(ndofn,force_rhs)
    end if

  end subroutine COU_RESIDUAL_FORCE

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Interpolation communication arrays
  !> @details Interpolation communication arrays
  !>          Input variables:
  !>          COUPLING ............................... Coupling structure
  !>          XX(NDIME,:) ............................ Coordinates of wet points
  !>          COUPLING % GEOME % NUMBER_WET_POINTS ... SIZE(XX,2)
  !>          COUPLING % ITYPE ....................... Vector projection
  !>          COUPLING % KIND ........................ BETWEEN_SUBDOMAINS/BETWEEN_ZONES
  !>          COUPLING % KDTREE ...................... Kdtree of source 
  !>          
  !>
  !----------------------------------------------------------------------
  
  subroutine COU_INIT_INTERPOLATE_POINTS_VALUES_NEW(coupling,COMM,CANDIDATE_SOURCE_NODES,CANDIDATE_SOURCE_ELEMENTS,CANDIDATE_SOURCE_BOUNDARIES)

    type(typ_color_coupling), intent(inout)          :: coupling
    integer(ip),              intent(in),   optional :: COMM
    logical(lg), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_NODES(:)
    integer(ip), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_ELEMENTS(:)
    logical(lg), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_BOUNDARIES(:)
    integer(ip)                                      :: icolo_source
    integer(ip)                                      :: icolo_target

    icolo_source = coupling % color_source
    icolo_target = coupling % color_target
    
    call COU_INIT_INTERPOLATE_POINTS_VALUES_OLD(coupling % wet % coord_wet,icolo_target,icolo_source,coupling,COMM,&
         CANDIDATE_SOURCE_NODES,CANDIDATE_SOURCE_ELEMENTS,CANDIDATE_SOURCE_BOUNDARIES)   
    
  end subroutine COU_INIT_INTERPOLATE_POINTS_VALUES_NEW
  

  subroutine COU_INIT_INTERPOLATE_POINTS_VALUES_OLD(xx,icolo,jcolo,coupling,COMM,CANDIDATE_SOURCE_NODES,CANDIDATE_SOURCE_ELEMENTS,CANDIDATE_SOURCE_BOUNDARIES)

    use mod_coupling_toolbox,  only : coupling_toolbox_points_in_partitions
    use mod_par_bin_structure, only : par_bin_structure
    use mod_par_bin_structure, only : par_bin_structure_deallocate
    use mod_par_bin_structure, only : par_bin_structure_initialization
    use mod_par_bin_structure, only : par_bin_structure_partition_bounding_box
    use def_coupli,            only : toler_relative_cou
    use def_coupli,            only : toler_absolute_cou
    use mod_parall,            only : typ_bin_structure
    use mod_parall,            only : par_color_to_code
    use mod_parall,            only : par_color_to_zone
    use mod_parall,            only : par_color_to_subd
    use mod_parall,            only : PAR_MY_WORLD_RANK 

    use mod_elmgeo,           only  : numb_natural_coordinates_bb      
    use mod_elmgeo,           only  : numb_natural_coordinates_ray      
    use mod_elmgeo,           only  : numb_natural_coordinates_nr     
    use mod_elmgeo,           only  : time_natural_coordinates_bb  
    use mod_elmgeo,           only  : time_natural_coordinates_ray 
    use mod_elmgeo,           only  : time_natural_coordinates_nr  
    use def_domain,           only  : nelty
    use def_domain,           only  : lexis
    use mod_htable,           only  : HtableMaxPrimeNumber
    use mod_htable,           only  : hash_t
    use mod_htable,           only  : htaini
    use mod_htable,           only  : htaadd
    use mod_htable,           only  : htades

   
    
    integer(ip),              intent(in)             :: icolo
    integer(ip),              intent(in)             :: jcolo
    real(rp),    pointer,     intent(in)             :: xx(:,:)
    type(typ_color_coupling), intent(inout)          :: coupling
    integer(ip),              intent(in),   optional :: COMM
    logical(lg), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_NODES(:)
    integer(ip), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_ELEMENTS(:)
    logical(lg), pointer,     intent(in),   optional :: CANDIDATE_SOURCE_BOUNDARIES(:)    
    integer(ip)                                      :: ii,jj,kk,pp,kpart,ipart,icoup
    integer(ip)                                      :: cpart_owner,number_wet_points
    integer(ip)                                      :: ielem,ierro,cz,ksize,kdtree_wet
    integer(ip)                                      :: cpart,ineig,isend,irecv,ielty
    integer(ip)                                      :: ksend,lsend,lrecv,krecv
    integer(ip)                                      :: kelem,inode,idime,inodb,pnodb
    integer(ip)                                      :: npoin_to_send,npoin_to_recv
    integer(ip)                                      :: nenti,kboun,iboun,ipoin,kpoin
    integer(ip)                                      :: jcoup,ipart_world
    integer(ip)                                      :: ipoin_min_max_send(2)
    integer(ip)                                      :: ipoin_min_max_recv(2)
    integer(ip)                                      :: ipoin_min,ipoin_max
    integer(ip)                                      :: PAR_CURRENT_SIZE,all_points_found
    integer(ip)                                      :: PAR_CURRENT_RANK,ipass
    integer(ip)                                      :: PAR_COMM_SAVE
    integer(4)                                       :: PAR_COMM_CURRENT4
    real(rp)                                         :: dista_min
    type(r2p),   pointer                             :: shapf(:)              
    integer(ip), pointer                             :: list_neighbors(:,:)
    logical(lg), pointer                             :: list_source_nodes(:)
    integer(ip), pointer                             :: npoin_send(:) 
    integer(ip), pointer                             :: npoin_recv(:) 
    integer(ip), pointer                             :: mask_npoin(:) 
    integer(ip), pointer                             :: check_points(:)
    logical(lg), pointer                             :: first_myself(:)
    type(i1p),   pointer                             :: decision_send(:) 
    type(i1p),   pointer                             :: decision_recv(:) 
    type(r1p),   pointer                             :: distance_send(:)      ! Distance check
    type(r1p),   pointer                             :: distance_recv(:)      ! Distance check
    real(rp),    pointer                             :: distance_min_value(:) ! Distance check
    integer(ip), pointer                             :: distance_min_cpart(:) ! Distance check
    type(i1p),   pointer                             :: my_part_to_point(:)
    type(i1p),   pointer                             :: my_point_to_part(:)
    type(r1p),   pointer                             :: coord_send(:)
    type(r1p),   pointer                             :: coord_recv(:)
    type(i1p),   pointer                             :: numer_send(:)
    type(i1p),   pointer                             :: numer_recv(:)
    logical(lg)                                      :: require_distance
    integer(ip), pointer                             :: lesou(:)             
    logical(lg), pointer                             :: lbsou(:)             
    logical(lg), pointer                             :: lnsou(:)             
    logical(lg), pointer                             :: intersection(:)      
    integer(ip), pointer                             :: PAR_WORLD_RANKS(:)   
    integer(4)                                       :: COMM4
    real(rp)                                         :: time0,time1,time2,time3,time4,time5
    type(hash_t)                                     :: ht
    !
    ! Bounding boxesec
    !
    type(typ_bin_structure)                          :: bin_structure
    real(rp)                                         :: comin(3),comax(3),delta(3)
    !
    ! For debugg
    !
    integer(ip)                                      :: lun_chec
    integer(ip)                                      :: lun_recv
    integer(ip), pointer                             :: num_chec(:)
    integer(ip), pointer                             :: num_chec_gat(:,:)
    integer(ip), pointer                             :: num_recv_gat(:,:)
    integer(4)                                       :: num_intersections4
    real(rp)                                         :: time10,time11
    integer(ip)                                      :: check_myself_first

    call cputim(time0)
    !
    ! Coupling
    !
    icoup        = coupling % number
    color_target = coupling % color_target
    color_source = coupling % color_source
    !
    ! Communicator
    !
    PAR_COMM_SAVE = PAR_COMM_CURRENT
    if( present(COMM) ) then
       PAR_COMM_CURRENT = COMM
    else
       PAR_COMM_CURRENT = PAR_COMM_COLOR(icolo,jcolo)
    end if
    PAR_COMM_CURRENT4 = int(PAR_COMM_CURRENT,4_4)
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_CURRENT,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    !
    ! Sizes 
    !
    if( .not. associated(xx) ) then
       number_wet_points = 0 
    else
       number_wet_points = size(xx,2)
       if( coupling % wet % number_wet_points /= 0 ) then
          if( coupling % wet % number_wet_points /= number_wet_points ) &
               call runend('SOMETHING STRANGE HAPPENS: '//intost(coupling % wet % number_wet_points)//', '//intost(number_wet_points))
       else
          coupling % wet % number_wet_points = number_wet_points 
       end if
    end if
    !
    ! Nullify and initialize
    !
    nullify(shapf)
    nullify(list_neighbors)
    nullify(npoin_send)
    nullify(npoin_recv)
    nullify(mask_npoin)
    nullify(check_points)
    nullify(first_myself)
    nullify(decision_send)
    nullify(decision_recv)
    nullify(distance_send)
    nullify(distance_recv)
    nullify(distance_min_value)
    nullify(distance_min_cpart)
    nullify(my_part_to_point)
    nullify(my_point_to_part)
    nullify(coord_send)
    nullify(coord_recv)
    nullify(numer_send)
    nullify(numer_recv)
    nullify(lesou)
    nullify(lbsou)
    nullify(lnsou)
    nullify(intersection)
    nullify(list_source_nodes)
    nullify(PAR_WORLD_RANKS)    
    comin = 0.0_rp
    comax = 0.0_rp
    cz    = PAR_CURRENT_SIZE
    !
    ! If the distance criterion is required. In the case of global numbering, this
    ! is uncessary, unless meshes are different!
    !
    if(     coupling % itype == NEAREST_BOUNDARY_NODE  .or. &
         &  coupling % itype == BOUNDARY_INTERPOLATION .or. &
         &  coupling % itype == ELEMENT_INTERPOLATION  .or. &
         &  coupling % itype == STRESS_PROJECTION      .or. &
         &  coupling % itype == PROJECTION             .or. &
         &  coupling % itype == NEAREST_ELEMENT_NODE   .or. &
         &  coupling % itype == BOUNDARY_VECTOR_PROJECTION ) then
       require_distance = .true.
    else
       require_distance = .false.
    end if
    !
    ! LESOU, LBSOU, LNSOU: Source elements/boundaries/nodes
    !
    if( INOTMASTER ) then

       if( present(CANDIDATE_SOURCE_ELEMENTS) ) then
          lesou => CANDIDATE_SOURCE_ELEMENTS
       else
          if(    coupling % itype == NEAREST_ELEMENT_NODE  .or. &
               & coupling % itype == NEAREST_BOUNDARY_NODE .or. &
               & coupling % itype == ELEMENT_INTERPOLATION ) then
             call memory_alloca(memor_cou,'NELEM','cou_init_interpolate_points_values',lesou,nelem)
             if( coupling % kind == BETWEEN_SUBDOMAINS ) then
                do ielem = 1,nelem
                   if( lesub(ielem) == coupling % subdomain_source ) then
                      lesou(ielem) = 1
                   end if
                end do
             else         
                lesou= 1
             end if
          end if
       end if

       if( present(CANDIDATE_SOURCE_BOUNDARIES) ) then
          lbsou => CANDIDATE_SOURCE_BOUNDARIES
       end if
          
       if( present(CANDIDATE_SOURCE_NODES) ) then
          lnsou => CANDIDATE_SOURCE_NODES
       else if( associated(lesou) ) then
          if(    coupling % itype == NEAREST_ELEMENT_NODE .or. &
               & coupling % itype == NEAREST_BOUNDARY_NODE ) then
             call memory_alloca(memor_cou,'NPOIN','cou_init_interpolate_points_values',lnsou,npoin)
             do ielem = 1,nelem
                if( lesou(ielem) == 1 ) then
                   do inode = 1,lnnod(ielem)
                      ipoin = lnods(inode,ielem)
                      lnsou(ipoin) = .true.
                   end do
                end if
             end do
          end if
       end if

    end if
    !
    ! Allocate memory
    !
    call memory_alloca(memor_cou,'NPOIN_SEND'      ,'cou_init_interpolate_points_values',npoin_send,      cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'NPOIN_RECV'      ,'cou_init_interpolate_points_values',npoin_recv,      cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'DECISION_SEND'   ,'cou_init_interpolate_points_values',decision_send,   cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'DECISION_RECV'   ,'cou_init_interpolate_points_values',decision_recv,   cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'MY_PART_TO_POINT','cou_init_interpolate_points_values',my_part_to_point,cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'MY_POINT_TO_PART','cou_init_interpolate_points_values',my_point_to_part,number_wet_points, 'INITIALIZE')
    call memory_alloca(memor_cou,'SHAPF'           ,'cou_init_interpolate_points_values',shapf,           cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'INTERSECTION'    ,'cou_init_interpolate_points_values',intersection,    cz,                'INITIALIZE',0_ip)
    call memory_alloca(memor_cou,'PAR_WORLD_RANKS' ,'cou_init_interpolate_points_values',PAR_WORLD_RANKS, cz,                'INITIALIZE',0_ip)

    if( coupling % itype == GLOBAL_NUMBERING ) then
       call memory_alloca(memor_cou,'NUMER_SEND','cou_init_interpolate_points_values',numer_send,      cz,                'INITIALIZE',0_ip)
       call memory_alloca(memor_cou,'NUMER_RECV','cou_init_interpolate_points_values',numer_recv,      cz,                'INITIALIZE',0_ip)       
    else
       call memory_alloca(memor_cou,'COORD_SEND','cou_init_interpolate_points_values',coord_send,      cz,                'INITIALIZE',0_ip)
       call memory_alloca(memor_cou,'COORD_RECV','cou_init_interpolate_points_values',coord_recv,      cz,                'INITIALIZE',0_ip)
    end if

    if( require_distance ) then
       call memory_alloca(memor_cou,'DISTANCE_SEND','cou_init_interpolate_points_values',distance_send,cz,'INITIALIZE',0_ip)
       call memory_alloca(memor_cou,'DISTANCE_RECV','cou_init_interpolate_points_values',distance_recv,cz,'INITIALIZE',0_ip)
    end if
    !
    ! World ranks. E.g. PAR_WORLD_RANKS(PAR_CURRENT_RANK) = PAR_MY_WORLD_RANK
    !
    if( present(COMM) ) then
       PAR_WORLD_RANKS(PAR_CURRENT_RANK) = PAR_MY_WORLD_RANK
       call PAR_MAX(PAR_CURRENT_SIZE,PAR_WORLD_RANKS(0:PAR_CURRENT_SIZE-1),PAR_COMM_IN4=PAR_COMM_CURRENT4)
       !cpart = PAR_CURRENT_RANK
       !PAR_WORLD_RANKS(cpart) = PAR_MY_WORLD_RANK
       !ipart = PAR_MY_WORLD_RANK
       !call PAR_ALLGATHER(ipart,PAR_WORLD_RANKS,PAR_COMM_IN=COMM)
       !print*,'www=',kfl_paral,PAR_WORLD_RANKS
    else
       do ipart = 0,PAR_WORLD_SIZE-1
          cpart = PAR_COMM_COLOR_PERM(icolo,jcolo,ipart)
          if( cpart > 0 ) PAR_WORLD_RANKS(cpart) = ipart
       end do
    end if
    !
    ! Bin structure for partitions
    !
    !call PAR_BARRIER() ; if(IMASTER) print*,'A'
    call par_bin_structure_initialization(bin_structure)
    call par_bin_structure_partition_bounding_box(&
         comin,comax,delta,&
         RELATIVE_TOLERANCE=toler_relative_cou,&
         ABSOLUTE_TOLERANCE=toler_absolute_cou)
    if( number_wet_points > 0 ) then
       do idime = 1,ndime
          comin(idime) = min(comin(idime),minval(xx(idime,1:number_wet_points)))
          comax(idime) = max(comax(idime),maxval(xx(idime,1:number_wet_points)))
       end do
    end if
    call par_bin_structure(bin_structure,comin,comax,COMM=PAR_COMM_CURRENT,VERBOSE=.false.)
    !
    ! Find list of partitions to which I must send my points
    ! PP:                            : 1 -> NUMBER_WET_POINTS
    ! NPOIN_SEND(CPART)              = number of wet points to send to CPART 
    ! MY_POINT_TO_PART(PP) % L(:)    = CPART, list of partitions where PP is in
    ! MY_PART_TO_POINT(CPART) % L(:) = PP,    list of points located in CPART bounding box
    !
    call coupling_toolbox_points_in_partitions(&
         number_wet_points,ndime,xx,jcolo,bin_structure,coupling % kfl_multi_source,&
         PAR_CURRENT_SIZE,PAR_CURRENT_RANK,PAR_WORLD_RANKS,&
         npoin_send,my_point_to_part,my_part_to_point)
    !
    ! Bounding box intersection. We will communicate only with these subdomains at first!
    !
    ! For vector projection, send the wet points to all the partitions as we have no way
    ! to know where the point falls. Idea: consider only the partitions bounding boxes
    ! crossed by the segment formed by the wet node ( x_wet , x_wet + 10000 * projection_vector )  
    ! 
    ! +-----+------+-------+-----+
    ! |     |      | x_wet |     |
    ! |     |      |    x  |     |
    ! +-----+------+----|--+-----+
    ! |     |      |    |  |     |
    ! |     |      |    |  |     |
    ! +-----+------+----|--+-----+
    ! |     |      |    |  |     |
    ! |     |      |    |  |     |
    ! +-----+------+----|--+-----+
    ! ///////////////// | ////////
    !                   |
    !                   | projection_vector
    !                   |
    !                   \/   
    !                   x 
    !
    do cpart = 0,PAR_CURRENT_SIZE-1
       intersection(cpart) = .true.
       do idime = 1,ndime
          if(    bin_structure % part_comin(idime,PAR_CURRENT_RANK) > bin_structure % part_comax(idime,cpart) .or. &
               & bin_structure % part_comin(idime,cpart)            > bin_structure % part_comax(idime,PAR_CURRENT_RANK) ) then
             intersection(cpart) = .false.
          end if
       end do
    end do
    if( coupling % itype == BOUNDARY_VECTOR_PROJECTION .and. INOTMASTER ) then
       do cpart = 1,PAR_CURRENT_SIZE-1
          intersection(cpart) = .true. 
       end do
    end if
    num_intersections4 = int(count(intersection,KIND=ip),KIND=4)
     
    !--------------------------------------------------------------------
    !
    !  1. How many wet points I send to (I) and how many I receive from (J)
    !
    !                            +-------+
    !                        =>  |       |  =>
    !                            |       |
    !          NPOIN_RECV(I) =>  |       |  => NPOIN_SEND(I)
    !                            +-------+
    !
    !  2. Send coordinates and receive coordinates
    !
    !                            +-------+
    !                        =>  |       |  =>
    !                            |       |
    !      COORD_RECV(I) % A =>  |       |  => COORD_SEND(I) % A
    !                            +-------+
    !
    !  3. Check if I have what I receive (find host element, nearest point, etc.) and send it back: DECISION_SEND(I) % L
    !     Receive the decision of the others: DECISION_RECV(I) % L
    !     In case of nearest point, send the minimal distance I have found as well: DISTANCE_SEND(I) % A
    !     Receive the minimal distance of others as well: DISTANCE_RECV(I) % A
    !
    !                            +-------+
    !                        <=  |       |  <=
    !                            |       |
    !   DECISION_SEND(I) % L <=  |       |  <= DECISION_RECV(I) % L
    !                            +-------+
    !                         
    !--------------------------------------------------------------------
    !
    ! Get how many points I should check if I have them
    !
    call PAR_ALLTOALL(1_ip,1_ip,npoin_send,npoin_recv,PAR_COMM_IN4=PAR_COMM_CURRENT4)
    call cputim(time1) 
    coupling % cputim(6) = coupling % cputim(6) + time1-time0
    !
    ! In the case of global numbering, reduce the possible number of communications
    ! In fact, the send and receive should be symmetric
    !
    if( coupling % itype == GLOBAL_NUMBERING .and. coupling % kfl_symmetry == 1 ) then
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( intersection(cpart) ) then
             if( npoin_send(cpart) == 0 .or. npoin_recv(cpart) == 0 ) then
                npoin_send(cpart)   = 0
                npoin_recv(cpart)   = 0
                intersection(cpart) = .false.
                call memory_deallo(memor_cou,'MY_PART_TO_POINT(CPART) % L', 'cou_init_interpolate_points_values',my_part_to_point(cpart) % l)
             elseif(1==2) then
                if( npoin_send(cpart) > 0 ) then
                   ipoin_min_max_send(1) = minval(lninv_loc(numer_send(cpart) % l))
                   ipoin_min_max_send(2) = maxval(lninv_loc(numer_send(cpart) % l))
                else
                   ipoin_min_max_send    = 0_rp
                end if
                call PAR_SEND_RECEIVE(2_ip,2_ip,ipoin_min_max_send,ipoin_min_max_recv,'IN CURRENT COUPLING',cpart,'BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                ipoin_min = max(ipoin_min_max_send(1),ipoin_min_max_recv(1))
                ipoin_max = min(ipoin_min_max_send(2),ipoin_min_max_recv(2))
                if( npoin_send(cpart) > 0_ip ) then
                   kk = 0
                   do ii = 1,npoin_to_send
                      pp    = my_part_to_point(cpart) % l(ii)
                      ipoin = coupling % wet % lpoin_wet(pp)
                      if( lninv_loc(ipoin) >= ipoin_min ) then
                         kk = kk + 1
                         my_part_to_point(cpart) % l(kk) = pp
                      end if
                   end do
                   if( kk /= ii ) then
                      npoin_send(cpart) = kk
                      call memory_resize(memor_cou,'MY_PART_TO_POINT(CPART) % L', 'cou_init_interpolate_points_values',my_part_to_point(cpart) % l,npoin_send(cpart))
                   end if
                end if
             end if
          end if
        end do
    end if
    !
    ! Save point coordinates
    !
    do cpart = 0,PAR_CURRENT_SIZE-1
       npoin_to_send = npoin_send(cpart)
       npoin_to_recv = npoin_recv(cpart)
       if( coupling % itype == GLOBAL_NUMBERING ) then
          call memory_alloca(memor_cou,'NUMER_SEND(CPART) % L', 'cou_init_interpolate_points_values',numer_send(cpart) % l,npoin_to_send)
          call memory_alloca(memor_cou,'NUMER_RECV(CPART) % L', 'cou_init_interpolate_points_values',numer_recv(cpart) % l,npoin_to_recv)
       else
          call memory_alloca(memor_cou,'COORD_SEND(CPART) % A', 'cou_init_interpolate_points_values',coord_send(cpart) % a,npoin_to_send*ndime)
          call memory_alloca(memor_cou,'COORD_RECV(CPART) % A', 'cou_init_interpolate_points_values',coord_recv(cpart) % a,npoin_to_recv*ndime)
       end if
       call memory_alloca(memor_cou,'DECISION_RECV(CPART) % L' ,'cou_init_interpolate_points_values',decision_recv(cpart)  % l,npoin_to_send)
       call memory_alloca(memor_cou,'DECISION_SEND(CPART) % L' ,'cou_init_interpolate_points_values',decision_send(cpart)  % l,npoin_to_recv)
       if( require_distance ) then
          call memory_alloca(memor_cou,'DISTANCE_SEND(CPART) % A','cou_init_interpolate_points_values',distance_send(cpart) % a,npoin_to_recv)
          call memory_alloca(memor_cou,'DISTANCE_RECV(CPART) % A','cou_init_interpolate_points_values',distance_recv(cpart) % a,npoin_to_send)
       end if
    end do
    !
    ! Copy points to send buffer
    !
    if( coupling % itype == GLOBAL_NUMBERING ) then
       do cpart = 0,PAR_CURRENT_SIZE-1
          npoin_to_send = npoin_send(cpart)
          if( npoin_to_send > 0 ) then
             do ii = 1,npoin_to_send
                pp    = my_part_to_point(cpart) % l(ii)
                ipoin = coupling % wet % lpoin_wet(pp)
                numer_send(cpart) % l(ii) = lninv_loc(ipoin)
             end do
          end if
       end do
    else
       do cpart = 0,PAR_CURRENT_SIZE-1
          npoin_to_send = npoin_send(cpart)
          if( npoin_to_send > 0 ) then
             kk = 0
             do ii = 1,npoin_to_send
                pp = my_part_to_point(cpart) % l(ii)
                do idime = 1,ndime
                   kk = kk + 1
                   coord_send(cpart) % a(kk) = xx(idime,pp)
                end do
             end do
          end if
       end do
    end if
    !
    ! Get points coordinates or node numbering
    !
    call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    if( coupling % itype == GLOBAL_NUMBERING ) then
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( intersection(cpart) ) then
             call PAR_SEND_RECEIVE(numer_send(cpart) % l,numer_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
          end if
       end do
    else
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( intersection(cpart) ) then
             call PAR_SEND_RECEIVE(coord_send(cpart) % a,coord_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
          end if
       end do
    end if

    call cputim(time2) 
    coupling % cputim(3) = coupling % cputim(3) + time2-time1
    call PAR_END_NON_BLOCKING_COMM(1_ip)
    call cputim(time2) 

    !---------------------------------------------------------------------------------------------
    !
    ! Look first for host elements in myslef if I am a source
    !
    !---------------------------------------------------------------------------------------------
        
    if( coupling % itype == ELEMENT_INTERPOLATION ) then
       
       check_myself_first = 0
       ipass              = 1
       
       if( par_part_in_color(PAR_MY_WORLD_RANK,jcolo) ) then
          !
          ! Check if I host the wet points
          ! CHECK_MYSLEF_FIRST= number of wet points I have found
          !       
          cpart = PAR_CURRENT_RANK
          if( npoin_recv(cpart) /= 0 ) then
             pp = npoin_recv(cpart)
             if( .not. associated(shapf(cpart)%a) ) &
                  call memory_alloca(memor_cou,'SHAPF(CPART) % A','cou_init_interpolate_points_values',shapf(cpart)%a,mnode,pp)
             call COU_WET_POINTS_HOST_ELEMENTS(ipass,pp,coord_recv(cpart) % a, distance_send(cpart) % a, shapf(cpart) % a,decision_send(cpart) % l,lesou)
             check_myself_first = count( decision_send(cpart) % l /= 0 , KIND=ip) 
          end if 
       end if
       
       call PAR_MAX(check_myself_first,PAR_COMM_IN4=PAR_COMM_CURRENT4)

       if( check_myself_first > 0 ) then
          !
          ! Mark my wet points and mark lost wet points not to be checked anymore by myself
          ! FIRST_MYSELF(PP) = .TRUE. if I found host element for PP
          !
          call memory_alloca(memor_cou,'FIRST_MYSELF','cou_init_interpolate_points_values',first_myself,number_wet_points)
          do ii = 1,npoin_recv(PAR_CURRENT_RANK)
             pp = my_part_to_point(PAR_CURRENT_RANK) % l(ii)
             if( decision_send(PAR_CURRENT_RANK) % l(ii) > 0 ) then
                first_myself(pp) = .true.                            ! I found PP!
             else
                decision_send(PAR_CURRENT_RANK) % l(ii) = my_huge    ! I did not find it, I will not check it anymore
             end if                
          end do
          !
          ! Tell people they do not have to look at the wet points I have already found
          !
          do cpart = 0,PAR_CURRENT_SIZE-1
             if( cpart /= PAR_CURRENT_RANK ) then
                do ii = 1,npoin_send(cpart)
                   pp = my_part_to_point(cpart) % l(ii)
                   if( first_myself(pp) ) decision_recv(cpart) % l(ii) = my_huge
                end do
             end if
          end do
          call memory_deallo(memor_cou,'FIRST_MYSELF','cou_init_interpolate_points_values',first_myself)
          !
          ! Send my first decision to my neighbor (put result in decision_send
          !
          call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
          call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
          do cpart = 0,PAR_CURRENT_SIZE-1
             if( intersection(cpart) .and. cpart /= PAR_CURRENT_RANK ) then
                call PAR_SEND_RECEIVE(decision_recv(cpart) % l,decision_send(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
             end if
          end do
          call cputim(time3) 
          coupling % cputim(4) = coupling % cputim(4) + time3-time2
          call PAR_END_NON_BLOCKING_COMM(1_ip)
          call cputim(time2) 
          
       end if
       
    end if
    
    !---------------------------------------------------------------------------------------------
    !
    ! Statistics, to be draw with gnuplot.
    ! **  is problem name
    ! *** is the coupling number
    !
    ! set xrange[0.5:]
    ! set yrange[0.5:]
    !
    ! set title  'Number of points rank x receives from rank y'
    ! set xlabel 'MPI rank x'
    ! set ylabel 'MPI rank y'
    !
    ! plot '**-histogram-recv-***.cou.res' matrix with image
    ! set title  'Number of points rank x receives from rank y and checked'
    ! set xlabel 'MPI rank x'
    ! set ylabel 'MPI rank y'
    ! plot '**-histogram-chec-***.cou.res' matrix with image
    !
    !---------------------------------------------------------------------------------------------
    
    if( 1 == 2 ) then
       nullify(num_chec)
       nullify(num_chec_gat)
       nullify(num_recv_gat)
       kpart = PAR_CURRENT_SIZE-1
       if( PAR_CURRENT_RANK == 0 ) then
          allocate(num_chec_gat(0:kpart,0:kpart))
          allocate(num_recv_gat(0:kpart,0:kpart))
       end if
       allocate(num_chec(0:kpart))
       num_chec = 0
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( npoin_recv(cpart) /= 0 ) then
             num_chec(cpart) = count(decision_send(cpart) % l == 0, KIND=ip)
          end if
       end do
       call PAR_GATHER(num_chec  ,num_chec_gat,'IN CURRENT COUPLING')
       call PAR_GATHER(npoin_recv,num_recv_gat,'IN CURRENT COUPLING') 
       if( PAR_CURRENT_RANK == 0 ) then
          lun_chec = iofile_available_unit()
          call iofile_open_unit(lun_chec,adjustl(trim(namda))//'-histogram-chec-'//trim(intost(icoup))//'.cou.res')
          lun_recv = iofile_available_unit()
          call iofile_open_unit(lun_recv,adjustl(trim(namda))//'-histogram-recv-'//trim(intost(icoup))//'.cou.res')
          do cpart = 0,kpart
             write(lun_chec,'(3000(1x,i6))') (num_chec_gat(cpart,ipart),ipart=0,kpart)
             write(lun_recv,'(3000(1x,i6))') (num_recv_gat(cpart,ipart),ipart=0,kpart)
             call iofile_flush_unit(lun_chec)
             call iofile_flush_unit(lun_recv)
          end do
          deallocate(num_chec_gat)
          deallocate(num_recv_gat)
          call iofile_close_unit(lun_chec)
          call iofile_close_unit(lun_recv)          
       end if
       deallocate(num_chec)
    end if
    
    !---------------------------------------------------------------------------------------------
    !
    ! Search strategy
    !
    !---------------------------------------------------------------------------------------------
    
    all_points_found = 0
    ipass            = 0
    call cputim(time2)
    if( coupling % itype == GLOBAL_NUMBERING ) then
       !
       ! Construct hash table
       !
       call htaini( ht, npoin, lidson=.true., AUTOMATIC_SIZE=.true.)
       if(associated(lninv_loc))then
          call htaadd( ht, npoin, lninv_loc)
       else
          if(npoin /= 0_ip) then
             call runend("Cou_init_interpolate_points_values: lninv_loc not associated")
          endif
       endif 
       if( ht % nelem /= npoin) call runend("Error: repited elements in lninv_loc") 
    endif
    
    do while( all_points_found == 0 )

       ipass = ipass + 1

       if( IMASTER    ) then

       else 

          if( coupling % itype == GLOBAL_NUMBERING .or. coupling % itype == SAME_COORDINATE ) then
             call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
          else
             call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
             call PAR_START_NON_BLOCKING_COMM(2_ip,num_intersections4)
          end if
          ! 
          ! Check if I own the points
          !       
          if(      coupling % itype == GLOBAL_NUMBERING ) then
             !
             ! Global numbering
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call COU_WET_POINTS_GLOBAL_NUMBERING(pp,numer_recv(cpart) % l, decision_send(cpart) % l, ht)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                end if
             end do

          else if( coupling % itype == ELEMENT_INTERPOLATION ) then
             !
             ! Element interpolation
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   if( .not. associated(shapf(cpart)%a) ) &
                        call memory_alloca(memor_cou,'SHAPF(CPART) % A','cou_init_interpolate_points_values',shapf(cpart)%a,mnode,pp)
                   call COU_WET_POINTS_HOST_ELEMENTS(ipass,pp,coord_recv(cpart) % a, distance_send(cpart) % a, shapf(cpart) % a,decision_send(cpart) % l,lesou)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(2_ip)
                   call PAR_SEND_RECEIVE(distance_send(cpart) % a,distance_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                end if
             end do

          else if( coupling % itype == NEAREST_BOUNDARY_NODE ) then
             !
             ! Nearest boundary node
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call COU_WET_POINTS_NEAREST_BOUNDARY_NODE(pp,coord_recv(cpart) % a,distance_send(cpart) % a,decision_send(cpart) % l,lnsou)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(2_ip)
                   call PAR_SEND_RECEIVE(distance_send(cpart) % a,distance_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                end if
             end do

          else if( coupling % itype == NEAREST_ELEMENT_NODE ) then
             !
             ! Nearest element node
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call COU_WET_POINTS_NEAREST_ELEMENT_NODE(pp,coord_recv(cpart) % a,distance_send(cpart) % a,decision_send(cpart) % l,lnsou)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(2_ip)
                   call PAR_SEND_RECEIVE(distance_send(cpart) % a,distance_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                end if
             end do

          else if( coupling % itype == SAME_COORDINATE ) then
             !
             ! Same coordinate
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call COU_WET_POINTS_SAME_COORDINATE(pp,coord_recv(cpart) % a,decision_send(cpart) % l,lnsou)
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                end if
             end do

          else if( coupling % itype == BOUNDARY_INTERPOLATION .or. &
               &   coupling % itype == STRESS_PROJECTION      .or. &
               &   coupling % itype == PROJECTION             ) then
             !
             ! Boundary interpolation and Gauss point interpolation
             ! Wet points are nodes in the first case and Gauss points in the second case
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call memory_alloca(memor_cou,'SHAPF(CPART) % A','cou_init_interpolate_points_values',shapf(cpart)%a,mnodb,pp)
                   call COU_WET_POINTS_HOST_BOUNDARIES(pp,coord_recv(cpart) % a,distance_send(cpart) % a,&
                        &                              decision_send(cpart) % l,shapf(cpart) % a,        &
                        &                              coupling % geome % kdtree,lbsou                   )
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(2_ip)
                   call PAR_SEND_RECEIVE(distance_send(cpart) % a,distance_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                end if
             end do

          else if( coupling % itype == BOUNDARY_VECTOR_PROJECTION ) then
             !
             ! Look for the boudnary corssing the projection along a vector
             !
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( npoin_recv(cpart) /= 0 ) then
                   pp = npoin_recv(cpart)
                   call memory_alloca(memor_cou,'SHAPF(CPART) % A','cou_init_interpolate_points_values',shapf(cpart)%a,mnodb,pp)
                   call COU_WET_POINTS_HOST_BOUNDARY_VECTOR(pp,coord_recv(cpart) % a,distance_send(cpart) % a,&
                        &                                   decision_send(cpart) % l,shapf(cpart) % a,        &
                        &                                   coupling % geome % kdtree                         )
                end if
                if( intersection(cpart) ) then
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
                   call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                   call PAR_SET_NON_BLOCKING_COMM_NUMBER(2_ip)
                   call PAR_SEND_RECEIVE(distance_send(cpart) % a,distance_recv(cpart) % a,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                end if
             end do

          end if
          if( coupling % itype == GLOBAL_NUMBERING .or. coupling % itype == SAME_COORDINATE ) then
             call PAR_END_NON_BLOCKING_COMM(1_ip)
          else
             call PAR_END_NON_BLOCKING_COMM(1_ip)
             call PAR_END_NON_BLOCKING_COMM(2_ip)
          end if
       end if
       !
       ! Check if we need a second round
       !
       if( coupling % itype == ELEMENT_INTERPOLATION .and. coupling % kfl_toda_costa == 1 ) then
          
          call memory_alloca(memor_cou,'CHECK_POINT','cou_init_interpolate_points_values',check_points,number_wet_points)
          !
          ! Check if someone has found the point
          !
          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart) 
                if( abs(decision_recv(cpart) % l(kk)) > 0 .and. decision_recv(cpart) % l(kk) < my_huge ) then
                   pp = my_part_to_point(cpart) % l(kk) 
                   check_points(pp) = 1
                end if
             end do
          end do
          if( number_wet_points > 0 ) then
             all_points_found = minval(check_points)
          else
             all_points_found = 1
          end if
          call PAR_MIN(all_points_found,'IN CURRENT COUPLING')
          !
          ! Tell neighbors which points they should look at 
          !
          !if(IMASTER) print*,'popo=',all_points_found
          if( all_points_found == 0 ) then
             !if( ipass == 2 ) call runend('WE ARE IN TROUBLE IN COUPLING')
             if( ipass == 2 ) all_points_found = 1
             call messages_live('WE GO FOR A SECOND ROUND TO FIND HOST ELEMENTS FOR SOME LOST WET NODES','WARNING')
             do cpart = 0,PAR_CURRENT_SIZE-1
                do kk = 1,npoin_send(cpart) 
                   pp = my_part_to_point(cpart) % l(kk)
                   if( check_points(pp) == 0 ) then
                      decision_recv(cpart) % l(kk) = -1
                   end if
                end do
             end do
             call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
             call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
             do cpart = 0,PAR_CURRENT_SIZE-1
                if( intersection(cpart) ) then
                   call PAR_SEND_RECEIVE(decision_recv(cpart) % l,decision_send(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
                end if
             end do
             call PAR_END_NON_BLOCKING_COMM(1_ip)
          end if

          call memory_deallo(memor_cou,'CHECK_POINTS','cou_init_interpolate_points_values',check_points)
       else
          all_points_found = 1
       end if

    end do
    if( coupling % itype == GLOBAL_NUMBERING ) then
       !
       ! Dellocate hast table
       ! 
       call htades( ht )
    
    endif
    
    !
    ! Remove -1
    !
    do cpart = 0,PAR_CURRENT_SIZE-1
       do kk = 1,npoin_send(cpart) 
          decision_recv(cpart) % l(kk) = max(decision_recv(cpart) % l(kk),0_ip)
          if( decision_recv(cpart) % l(kk) >= my_huge ) decision_recv(cpart) % l(kk) = 0
       end do
       do kk = 1,npoin_recv(cpart) 
          decision_send(cpart) % l(kk) = max(decision_send(cpart) % l(kk),0_ip)
          if( decision_send(cpart) % l(kk) >= my_huge ) decision_send(cpart) % l(kk) = 0
       end do
    end do


    if( INOTMASTER ) then
       !
       ! Decide the owner 
       !
       if( ( coupling % itype == GLOBAL_NUMBERING .or. coupling % itype == SAME_COORDINATE ) .and. number_wet_points > 0 .and. coupling % kfl_multi_source == 0 ) then
          !
          ! For element interpolation take element which belongs to subdomain
          ! with highest rank
          !
          call memory_alloca(memor_cou,'DISTANCE_MIN_CPART','cou_init_interpolate_points_values',distance_min_cpart,number_wet_points)

          distance_min_cpart = 0_ip

          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart) 
                if( decision_recv(cpart) % l(kk) /= 0 ) then
                   pp = my_part_to_point(cpart) % l(kk) 
                   if( cpart > distance_min_cpart(pp) ) then
                      distance_min_cpart(pp) = cpart
                   end if
                end if
             end do
          end do

          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart) 
                pp = my_part_to_point(cpart) % l(kk)
                if( distance_min_cpart(pp) /= cpart ) decision_recv(cpart) % l(kk) = 0
             end do
          end do

          do pp = 1,number_wet_points
             if( associated(my_point_to_part(pp) % l) ) my_point_to_part(pp) % l(1) = distance_min_cpart(pp)
          end do

          call memory_deallo(memor_cou,'DISTANCE_MIN_CPART','cou_init_interpolate_points_values',distance_min_cpart)

       else if( require_distance .and. number_wet_points > 0 ) then
          !
          ! Look for minimum distance
          ! 1. Compute minimum distance for each wet point pp, from 1 to number_wet_points
          ! 2. Put decision_recv(cpart) % l(kk) = 0 when cpart is different from the owner
          ! 3. my_point_to_part(pp) % l(1) = cpart_owner
          ! 
          call memory_alloca(memor_cou,'DISTANCE_MIN_VALUE','cou_init_interpolate_points_values',distance_min_value,number_wet_points)
          call memory_alloca(memor_cou,'DISTANCE_MIN_CPART','cou_init_interpolate_points_values',distance_min_cpart,number_wet_points)

          distance_min_value = huge(1.0_rp)

          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart) 
                if( decision_recv(cpart) % l(kk) /= 0 ) then
                   pp = my_part_to_point(cpart) % l(kk) 
                   if( distance_recv(cpart) % a(kk) < distance_min_value(pp) ) then
                      distance_min_value(pp)      = distance_recv(cpart) % a(kk)
                      distance_min_cpart(pp)      = cpart
                   end if
                end if
             end do
          end do

          do cpart = 0,PAR_CURRENT_SIZE-1
             do kk = 1,npoin_send(cpart) 
                pp = my_part_to_point(cpart) % l(kk)
                if( distance_min_cpart(pp) /= cpart ) decision_recv(cpart) % l(kk) = 0
             end do
          end do

          do pp = 1,number_wet_points
             if( associated(my_point_to_part(pp) % l) ) my_point_to_part(pp) % l(1) = distance_min_cpart(pp)
          end do

          call memory_deallo(memor_cou,'DISTANCE_MIN_CPART','cou_init_interpolate_points_values',distance_min_cpart)
          call memory_deallo(memor_cou,'DISTANCE_MIN_VALUE','cou_init_interpolate_points_values',distance_min_value)

       end if
       !
       ! Send my result back
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,num_intersections4)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do cpart = 0,PAR_CURRENT_SIZE-1
          if( intersection(cpart) ) then
             call PAR_SEND_RECEIVE(decision_recv(cpart) % l,decision_send(cpart) % l,'IN CURRENT COUPLING',cpart,'NON BLOCKING',PAR_COMM_IN4=PAR_COMM_CURRENT4)
          end if
       end do
       call cputim(time3) 
       coupling % cputim(4) = coupling % cputim(4) + time3-time2
       call PAR_END_NON_BLOCKING_COMM(1_ip)
       call cputim(time3) 

    end if

    if( INOTMASTER ) then

       !-----------------------------------------------------------------
       !
       ! decision_send(cpart) % l(1:npoin_recv(cpart)) /= 0 => I am in charge of this point OF CPART
       ! decision_recv(cpart) % l(1:npoin_send(cpart)) /= 0 => I need point from CPART
       !
       !-----------------------------------------------------------------
       !
       ! Fill in type
       !
       call memory_alloca(memor_cou,'LIST_NEIGHBORS','cou_init_interpolate_points_values',list_neighbors,2_ip,PAR_CURRENT_SIZE,'INITIALIZE',1_ip,0_ip)
       !
       ! Number of subdomains to communicate with
       !
       coupling % commd % nneig = 0
       do cpart = 0,PAR_CURRENT_SIZE-1
          do ii = 1,npoin_recv(cpart)
             list_neighbors(1,cpart) = list_neighbors(1,cpart) + min(1_ip,decision_send(cpart) % l(ii)) 
          end do
          do ii = 1,npoin_send(cpart)
             list_neighbors(2,cpart) = list_neighbors(2,cpart) + min(1_ip,decision_recv(cpart) % l(ii))
          end do
          coupling % commd % nneig = coupling % commd % nneig &
               + min(max(list_neighbors(1,cpart),list_neighbors(2,cpart)),1_ip)
       end do
       !
       ! Send and receives sizes of subdomains
       ! COMMD % LSEND_SIZE(:)
       ! COMMD % LRECV_SIZE(:)
       !
       call memory_alloca(memor_cou,'COMMD % NEIGHTS',   'cou_init_interpolate_points_values',coupling % commd % neights,   coupling % commd % nneig)
       call memory_alloca(memor_cou,'COMMD % LSEND_SIZE','cou_init_interpolate_points_values',coupling % commd % lsend_size,coupling % commd % nneig+1)
       call memory_alloca(memor_cou,'COMMD % LRECV_SIZE','cou_init_interpolate_points_values',coupling % commd % lrecv_size,coupling % commd % nneig+1)
       ineig = 0
       do cpart = 0,PAR_CURRENT_SIZE-1
          isend = list_neighbors(1,cpart)
          irecv = list_neighbors(2,cpart)
          if( max(isend,irecv) > 0 ) then
             ineig = ineig + 1
             coupling % commd % neights(ineig)    = cpart
             coupling % commd % lsend_size(ineig) = isend
             coupling % commd % lrecv_size(ineig) = irecv
          end if
       end do
       !
       ! Size to list
       !
       ksend = coupling % commd % lsend_size(1)
       krecv = coupling % commd % lrecv_size(1)
       coupling % commd % lsend_size(1) = 1
       coupling % commd % lrecv_size(1) = 1
       do ineig = 2,coupling % commd % nneig + 1
          lsend = coupling % commd % lsend_size(ineig)
          lrecv = coupling % commd % lrecv_size(ineig)
          coupling % commd % lsend_size(ineig) = coupling % commd % lsend_size(ineig-1) + ksend
          coupling % commd % lrecv_size(ineig) = coupling % commd % lrecv_size(ineig-1) + krecv
          ksend = lsend
          krecv = lrecv
       end do
       coupling % commd % lsend_dim = coupling % commd % lsend_size(coupling % commd % nneig + 1) - 1
       coupling % commd % lrecv_dim = coupling % commd % lrecv_size(coupling % commd % nneig + 1) - 1
       !
       ! Order points
       ! KK is the order I receive 
       ! 
       call memory_alloca(memor_cou,'STATUS','cou_init_interpolate_points_values',coupling % geome % status,coupling % wet % number_wet_points)                    
       kk = 0

       do cpart = 0,PAR_CURRENT_SIZE-1
          do ii = 1,npoin_send(cpart)
             pp = my_part_to_point(cpart) % l(ii)
             if( associated(my_point_to_part(pp) % l) ) then
                if( my_point_to_part(pp) % l(1) == cpart ) then
                   !
                   ! CPART has wet point pp
                   !
                   kk = kk + 1
                   coupling % geome % status(pp) = kk
                end if
             end if
          end do
       end do
       !
       ! Allocate geometrical information
       ! NENTI= number of entities (elements/boundaries/nodes)
       !
       nenti = 0
       do cpart = 0,PAR_CURRENT_SIZE-1
          nenti = nenti + list_neighbors(1,cpart)
       end do

       if( coupling % itype == ELEMENT_INTERPOLATION ) then
          !
          ! Element interpolation
          !
          coupling % geome % nelem_source = nenti
          call memory_alloca(memor_cou,'GEOME % SHAPF'            ,'cou_init_interpolate_points_values',coupling % geome % shapf,mnode, coupling % geome % nelem_source)
          call memory_alloca(memor_cou,'GEOME % LELEM'            ,'cou_init_interpolate_points_values',coupling % geome % lelem_source,coupling % geome % nelem_source)
          call memory_alloca(memor_cou,'GEOME % LIST_SOURCE_NODES','cou_init_interpolate_points_values',list_source_nodes,npoin)

          kelem = 0
          do cpart = 0,PAR_CURRENT_SIZE-1
             kk = 0
             do ii = 1,npoin_recv(cpart)
                ielem = decision_send(cpart) % l(ii)
                if( ielem > 0 ) then
                   kelem = kelem + 1
                   coupling % geome % lelem_source(kelem)  = ielem
                   do inode = 1,mnode
                      coupling % geome % shapf(inode,kelem) = shapf(cpart) % a(inode,ii)
                   end do
                   do inode = 1,lnnod(ielem)
                      ipoin = lnods(inode,ielem)
                      list_source_nodes(ipoin) = .true.
                   end do
                end if
             end do
          end do

       else if( coupling % itype == NEAREST_BOUNDARY_NODE      .or. &
            &   coupling % itype == NEAREST_ELEMENT_NODE       .or. &
            &   coupling % itype == BOUNDARY_VECTOR_PROJECTION .or. &
            &   coupling % itype == GLOBAL_NUMBERING           .or. &
            &   coupling % itype == SAME_COORDINATE            ) then
          !
          ! Nearest boundary node/global numbering
          !
          coupling % geome % npoin_source = nenti          
          call memory_alloca(memor_cou,'GEOME % LPOIN_SOURCE','cou_init_interpolate_points_values',coupling % geome % lpoin_source,coupling % geome % npoin_source)  

          kpoin = 0
          do cpart = 0,PAR_CURRENT_SIZE-1
             kk = 0
             do ii = 1,npoin_recv(cpart)
                ipoin = decision_send(cpart) % l(ii)
                if( ipoin > 0 ) then
                   kpoin = kpoin + 1
                   coupling % geome % lpoin_source(kpoin)  = ipoin
                end if
             end do
          end do

       else if( coupling % itype == BOUNDARY_INTERPOLATION ) then       
          !
          ! Boundary interpolation
          !  
          coupling % geome % nboun_source = nenti          
          call memory_alloca(memor_cou,'GEOME % SHAPF'            ,'cou_init_interpolate_points_values',coupling % geome % shapf,mnodb, coupling % geome % nboun_source)
          call memory_alloca(memor_cou,'GEOME % LBOUN'            ,'cou_init_interpolate_points_values',coupling % geome % lboun_source,coupling % geome % nboun_source)
          call memory_alloca(memor_cou,'LIST_SOURCE_NODES','cou_init_interpolate_points_values',list_source_nodes,npoin)

          kboun = 0
          do cpart = 0,PAR_CURRENT_SIZE-1
             kk = 0
             do ii = 1,npoin_recv(cpart)
                iboun = decision_send(cpart) % l(ii)
                if( iboun > 0 ) then
                   kboun = kboun + 1
                   coupling % geome % lboun_source(kboun)  = iboun
                   do inodb = 1,mnodb
                      coupling % geome % shapf(inodb,kboun) = shapf(cpart) % a(inodb,ii)
                   end do
                   do inodb = 1,lnnob_cou(iboun)
                      ipoin = lnodb_cou(inodb,iboun)
                      list_source_nodes(ipoin) = .true.
                   end do
                end if
             end do
          end do

       else if( coupling % itype == STRESS_PROJECTION .or. coupling % itype == PROJECTION ) then       
          !
          ! Projection
          !  
          coupling % geome % nboun_source = nenti          
          call memory_alloca(memor_cou,'GEOME % SHAPF',            'cou_init_interpolate_points_values',coupling % geome % shapf,mnodb,coupling % geome % nboun_source)
          call memory_alloca(memor_cou,'GEOME % LBOUN',            'cou_init_interpolate_points_values',coupling % geome % lboun_source,coupling % geome % nboun_source)
          call memory_alloca(memor_cou,'LIST_SOURCE_NODES','cou_init_interpolate_points_values',list_source_nodes,npoin)

          kboun = 0
          do cpart = 0,PAR_CURRENT_SIZE-1
             kk = 0
             do ii = 1,npoin_recv(cpart)
                iboun = decision_send(cpart) % l(ii)
                if( iboun > 0 ) then
                   kboun = kboun + 1
                   coupling % geome % lboun_source(kboun)  = iboun
                   do inodb = 1,mnodb
                      coupling % geome % shapf(inodb,kboun) = shapf(cpart) % a(inodb,ii)
                   end do
                   do inodb = 1,lnnob_cou(iboun)
                      ipoin = lnodb_cou(inodb,iboun)
                      list_source_nodes(ipoin) = .true.
                   end do
                   !write(*,'(i2,2(1x,e12.6),i5,i5)') kfl_paral,coord_recv(cpart) % a((ii-1)*2+1),coord_recv(cpart) % a((ii-1)*2+2),lninv_loc(lnodb(1:2,iboun))
                end if
             end do
          end do
       end if
       !
       ! Allocate sched_perm for scheduling permutation
       !
       call memory_alloca(memor_cou,'GEOME % SCHED_PERM','cou_init_interpolate_points_values',coupling % geome % sched_perm,coupling % commd % lrecv_dim)
       !
       ! Count number of source nodes if not already done
       !
       if( associated(list_source_nodes) ) then
          coupling % geome % npoin_source = 0
          do ipoin = 1,npoin
             if( list_source_nodes(ipoin) ) coupling % geome % npoin_source = coupling % geome % npoin_source + 1
          end do
          if( associated(coupling % geome % lpoin_source) ) then
             call runend('POINTER ALREADY ASSOCIATED')
          else
             call memory_alloca(memor_cou,'GEOME % LPOIN','cou_init_interpolate_points_values',coupling % geome % lpoin_source,coupling % geome % npoin_source)  
          end if
          kpoin = 0
          do ipoin = 1,npoin
             if( list_source_nodes(ipoin) ) then
                kpoin = kpoin + 1
                coupling % geome % lpoin_source(kpoin) = ipoin
             end if
          end do
       end if
       !
       ! Allocate memory for previous values
       !
       !call memory_alloca(memor_cou,'VALUES','cou_init_interpolate_points_values',coupling % values,coupling % geome % number_wet_points)
       !
       ! Check points that do not have a host partition
       ! If receive:        COUPLI % GEOME % LRECV_DIM values from other subdomains
       ! If need values on: COUPLI % GEOME % NUMBER_WET_POINTS values
       ! COUPLI % GEOME % STATUS(IPOIN) = JPOIN (from 1 to COUPLI % GEOME % LRECV_DIM)
       !
       ierro = 0
       ii    = 0
       do pp = 1,number_wet_points
          if( associated(my_point_to_part(pp) % l) ) then
             if( my_point_to_part(pp) % l(1) == 0 ) then
                ierro = 1
             else
                ii = ii + 1
             end if
          else
             ierro = 1
          end if
       end do

    else
       call cputim(time3)
    end if 
    !
    ! Allocate coupling % values_frequ for FSI in case of frequency is active for exchange
    !
    if( coupling % frequ_send > 1_ip .or. coupling % frequ_recv > 1_ip ) then
       call memory_alloca(memor_cou,'VALUES_FREQU','mod_couplings',coupling % values_frequ,ndime,coupling % geome % npoin_source,2_ip)
       do ipoin = 1_ip, coupling % geome % npoin_source
          do idime = 1_ip, ndime
             coupling % values_frequ(idime, ipoin, :) = 0_rp
          end do
       end do
    end if
    
    call cputim(time4)
    coupling % cputim(5) = coupling % cputim(5) + time4-time3
    !
    ! Deallocate memory
    !
    call memory_deallo(memor_cou,'LIST_NEIGHBORS'    ,'cou_init_interpolate_points_values',list_neighbors)
    call memory_deallo(memor_cou,'LIST_SOURCE_NODES' ,'cou_init_interpolate_points_values',list_source_nodes)
    call memory_deallo(memor_cou,'NPOIN_SEND'        ,'cou_init_interpolate_points_values',npoin_send)
    call memory_deallo(memor_cou,'NPOIN_RECV'        ,'cou_init_interpolate_points_values',npoin_recv)
    call memory_deallo(memor_cou,'MASK_NPOIN'        ,'cou_init_interpolate_points_values',mask_npoin)
    call memory_deallo(memor_cou,'DECISION_SEND'     ,'cou_init_interpolate_points_values',decision_send)
    call memory_deallo(memor_cou,'DECISION_RECV'     ,'cou_init_interpolate_points_values',decision_recv)
    call memory_deallo(memor_cou,'COORD_SEND'        ,'cou_init_interpolate_points_values',coord_send)
    call memory_deallo(memor_cou,'COORD_RECV'        ,'cou_init_interpolate_points_values',coord_recv)   
    call memory_deallo(memor_cou,'MY_PART_TO_POINT'  ,'cou_init_interpolate_points_values',my_part_to_point)
    call memory_deallo(memor_cou,'MY_POINT_TO_PART'  ,'cou_init_interpolate_points_values',my_point_to_part)
    call memory_deallo(memor_cou,'SHAPF'             ,'cou_init_interpolate_points_values',shapf)
    call memory_deallo(memor_cou,'DISTANCE_SEND'     ,'cou_init_interpolate_points_values',distance_send)
    call memory_deallo(memor_cou,'DISTANCE_RECV'     ,'cou_init_interpolate_points_values',distance_recv)
    call memory_deallo(memor_cou,'INTERSECTION'      ,'cou_init_interpolate_points_values',intersection)
    call memory_deallo(memor_cou,'PAR_WORLD_RANKS'   ,'cou_init_interpolate_points_values',PAR_WORLD_RANKS)
    call memory_deallo(memor_cou,'NUMER_SEND'        ,'cou_init_interpolate_points_values',numer_send)
    call memory_deallo(memor_cou,'NUMER_RECV'        ,'cou_init_interpolate_points_values',numer_recv)
    if( .not. present(CANDIDATE_SOURCE_ELEMENTS) ) then
       call memory_deallo(memor_cou,'LESOU     '     ,'cou_init_interpolate_points_values',lesou)
    end if
    if( .not. present(CANDIDATE_SOURCE_NODES) ) then
       call memory_deallo(memor_cou,'LNSOU     '     ,'cou_init_interpolate_points_values',lnsou)
    end if
   if( .not. present(CANDIDATE_SOURCE_BOUNDARIES) ) then
       call memory_deallo(memor_cou,'LBSOU     '     ,'cou_init_interpolate_points_values',lbsou)
    end if

    call par_bin_structure_deallocate(bin_structure)
    !
    ! Recover communicator 
    !
    PAR_COMM_CURRENT = PAR_COMM_SAVE

    !if( PAR_MY_WORLD_RANK == 0 ) then
    !   do ii=1,nx;do jj=1,ny ; do kk=1,nz
    !      if( bin_size(ii,jj,kk) > 0 ) then
    !         do iboxe = 1,bin_size(ii,jj,kk)
    !            ipart = par_bin_part(ii,jj,kk) % l(iboxe)
    !            icode = PAR_COMM_WORLD_TO_CODE_PERM(1,ipart)
    !            ipart = PAR_COMM_WORLD_TO_CODE_PERM(2,ipart)
    !         end do
    !      end if
    !   end do; end do; end do
    !end if
    !
    ! Timings
    !
    if( 1 == 2 ) then
       do ielty = 1,nelty
          if( lexis(ielty) == 1 ) then
             write(kfl_paral+90,*) '---------------------------------'
             write(kfl_paral+90,*) 'TYPE=',ielty
             write(kfl_paral+90,*) numb_natural_coordinates_bb(ielty)     
             write(kfl_paral+90,*) numb_natural_coordinates_ray(ielty)           
             write(kfl_paral+90,*) numb_natural_coordinates_nr(ielty)          
             write(kfl_paral+90,*) time_natural_coordinates_bb(ielty)       
             write(kfl_paral+90,*) time_natural_coordinates_ray(ielty)      
             write(kfl_paral+90,*) time_natural_coordinates_nr(ielty)
             call iofile_flush_unit(kfl_paral+90_ip)
          end if
       end do
    end if

    !call PAR_MAX(20_ip,coupling % cputim)
    !if(IMASTER) print*,coupling % cputim(1:8)
    !call runend('O.K.!')
  end subroutine COU_INIT_INTERPOLATE_POINTS_VALUES_OLD

  subroutine COU_PRESCRIBE_DIRICHLET_IN_MATRIX(nbvar,npopo,ia,ja,an,bb,xx)
    integer(ip), intent(in)    :: nbvar
    integer(ip), intent(in)    :: npopo
    integer(ip), intent(in)    :: ia(*)
    integer(ip), intent(in)    :: ja(*)
    real(rp),    intent(inout) :: an(nbvar,nbvar,*)
    real(rp),    intent(inout) :: xx(nbvar,*)
    real(rp),    intent(inout) :: bb(nbvar,*)
    integer(ip)                :: nrows,icoup,ii,jj,kk,nn,ll
    integer(ip)                :: izdom
    real(rp),    pointer       :: xx_tmp(:,:)

    if(    INOTMASTER .and. ( &
         & I_AM_INVOLVED_IN_A_COUPLING_TYPE(BETWEEN_SUBDOMAINS,DIRICHLET_IMPLICIT) .or. &
         & I_AM_INVOLVED_IN_A_COUPLING_TYPE(BETWEEN_SUBDOMAINS,DIRICHLET_EXPLICIT) ) ) then

       nrows = nbvar * npopo
       if( nbvar > 1 ) call runend('COU_PRESCRIBE_DIRICHLET_IN_MATRIX: NOT CODED')
       nullify( xx_tmp )
       allocate( xx_tmp(nbvar,npopo) )
       do ii = 1,npopo
          do nn = 1,nbvar
             xx_tmp(nn,ii) = xx(nn,ii)
          end do
       end do
       do icoup = 1,mcoup
          if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
             if( coupling_type(icoup) % what == DIRICHLET_EXPLICIT .or. coupling_type(icoup) % what == DIRICHLET_IMPLICIT ) then
                call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,xx,xx_tmp)
             end if
          end if
       end do
       deallocate( xx_tmp )
       do icoup = 1,mcoup
          if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
             if( coupling_type(icoup) % what == DIRICHLET_EXPLICIT .or. coupling_type(icoup) % what == DIRICHLET_IMPLICIT ) then
                do kk = 1,coupling_type(icoup) % wet % npoin_wet
                   ii  = coupling_type(icoup) % wet % lpoin_wet(kk)
                   do izdom = ia(ii),ia(ii+1)-1
                      jj = ja(izdom)
                      if( ii == jj ) then
                         an(1,1,izdom) = 1.0_rp
                      else
                         an(1,1,izdom) = 0.0_rp                    
                      end if
                   end do
                   bb(1,ii) = xx(1,ii)
                end do
             end if
          end if
       end do
    end if

  end subroutine COU_PRESCRIBE_DIRICHLET_IN_MATRIX

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check who owns the same node with same global numbering
  !> @details For each partition whcih have received some points,
  !>          check if I have a node with same global numbering
  !
  !----------------------------------------------------------------------
  subroutine COU_WET_POINTS_GLOBAL_NUMBERING(npoin_test,numer_test,inout_test, ht)

    use mod_htable

    implicit none
    integer(ip),          intent(in)  :: npoin_test                    !< Number of points to test
    integer(ip),          intent(in)  :: numer_test(npoin_test)        !< Global numbering of test points
    integer(ip),          intent(out) :: inout_test(npoin_test)        !< In or outside en element
     type(hash_t),          intent(in)  :: ht

    integer(ip)                       :: kpoin,ipoin_local
    integer(ip)                       :: ipoin_global
    
    !
    ! obtain lids
    !
    inout_test(:) = 0_ip
     
     do kpoin = 1,npoin_test
        ipoin_global = numer_test(kpoin)
       ipoin_local  = htalid( ht,ipoin_global)
        if( ipoin_global > 0 ) then
           inout_test(kpoin) = ipoin_local
       endif
    end do

     ! do kpoin = 1,npoin
     !    ipoin_global = lninv_loc(kpoin)
     !    ipoin_local  = htalid( ht,ipoin_global)
     !    if( ipoin_local > 0 ) then
     !       inout_test(ipoin_local) = kpoin
     !    endif
     ! end do


  
 end subroutine COU_WET_POINTS_GLOBAL_NUMBERING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          check which of them host them
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_HOST_ELEMENTS(ipass,npoin_test,coord_test,dimin_test,shapf_test,inout_test,lesou)

    use def_domain, only : lelbo
    use def_master, only : kfl_paral
    use mod_elsest
    integer(ip),          intent(in)   :: npoin_test                    !< Number of points to test
    real(rp),             intent(in)   :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    real(rp),             intent(out)  :: dimin_test(npoin_test)        !< Minimum distance to a boundary node
    real(rp),             intent(out)  :: shapf_test(mnode,npoin_test)  !< Shape function in elements
    integer(ip),          intent(out)  :: inout_test(npoin_test)        !< In or outside en element
    integer(ip), pointer, intent(in)   :: lesou(:)
    integer(ip),          intent(in)   :: ipass
    integer(ip)                        :: kpoin,ielem,dummi(2)
    integer(ip)                        :: iboun,pnodb,ii
    integer(ip)                        :: ielse15,ielse8,ielse14
    real(rp)                           :: coloc(3),deriv(64)
    real(rp)                           :: dista
    integer(ip), pointer               :: my_lesou(:)
    real(rp)                           :: time1,time2
    
    nullify(my_lesou)
    !
    ! Tell elsest that not all the elements should be considered
    ! Consider element only if LESOU(IELEM) = 1
    !
    ielse14 = ielse(14)  
    ielse15 = ielse(15)  
    ielse8  = ielse( 8)
    !
    ! Second chance, use A_TODA_COSTA options, but looking at nearest
    ! boundary elements
    !
    if( ipass == 2 ) then
       ielse(8) = ELSEST_A_TODA_COSTA
       call memory_alloca(memor_cou,'MY_LESOU','COU_WET_POINTS_HOST_ELEMENTS',my_lesou,nelem,'INITIALIZE')
       do iboun = 1,nboun
          pnodb = lnnob(iboun)
          ielem = lelbo(iboun)
          my_lesou(ielem) = 1
       end do
       do ielem = 1,nelem
          my_lesou(ielem) = min(my_lesou(ielem),lesou(ielem))
       end do
       ielse(14) = 1
       ielse(15) = 1
    else
       my_lesou => lesou
    end if
    !
    ! Loop over test points
    !
    !---------------------------------------------------------------------------     
    !$OMP PARALLEL DO                                                          &
    !$OMP SCHEDULE     ( DYNAMIC , 100 )                                       & 
    !$OMP DEFAULT      ( NONE )                                                &
    !$OMP SHARED       ( dimin_test, ielse, relse, meshe, coord_test,          &
    !$OMP                shapf_test, my_lesou, inout_test, ndivi, npoin_test ) &    
    !$OMP PRIVATE      ( kpoin, ielem, deriv, coloc, dista )                                            
    !---------------------------------------------------------------------------          
    
    do kpoin = 1,npoin_test
       
       if( inout_test(kpoin) <= 0 ) then
          
          ielem = 0
          dimin_test(kpoin) = huge(1.0_rp)

          call elsest_host_element(&
               ielse,relse,1_ip,meshe(ndivi),coord_test(:,kpoin),ielem,&
               shapf_test(:,kpoin),deriv,coloc,dista,my_lesou)

          if( ielem <= 0 ) then
             inout_test(kpoin) = 0              
          else
             if( my_lesou(ielem) == 1 ) then
                inout_test(kpoin) = ielem
                dimin_test(kpoin) = dista
             else
                inout_test(kpoin) = 0
             end if
          end if

       else if( inout_test(kpoin) >= huge(1_ip) ) then

          inout_test(kpoin) = my_huge
          dimin_test(kpoin) = huge(1.0_rp)
                
       end if
       
    end do
    !$OMP END PARALLEL DO 

    ielse( 8) = ielse8
    ielse(14) = ielse14
    ielse(15) = ielse15
    
    if( ipass == 2 ) then
       ielse(8) = ielse8
       call memory_deallo(memor_cou,'MY_LESOU','COU_WET_POINTS_HOST_ELEMENTS',my_lesou)
    end if

  end subroutine COU_WET_POINTS_HOST_ELEMENTS

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          check which of them host them
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_NEAREST_BOUNDARY_NODE(npoin_test,coord_test,dimin_test,inout_test,lnsou)
    integer(ip), intent(in)          :: npoin_test                    !< Number of points to test
    real(rp),    intent(in)          :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    real(rp),    intent(out)         :: dimin_test(npoin_test)        !< Minimum distance to a boundary node
    integer(ip), intent(out)         :: inout_test(npoin_test)        !< Node with minimum distance
    logical(lg), intent(in), pointer :: lnsou(:)
    integer(ip)                      :: ipoin,ibono,ipoin_min
    integer(ip)                      :: idime,kpoin
    real(rp)                         :: dista,dista_min

    do kpoin = 1,npoin_test

       dista_min = huge(1.0_rp)
       ipoin_min = 0

       do ibono = 1,nbono
          ipoin = lbono(ibono)
          if( lnsou(ipoin) ) then
             dista = 0.0_rp
             do idime = 1,ndime
                dista = dista + ( coord_test(idime,kpoin) - coord(idime,ipoin) ) ** 2
             end do
             if( dista < dista_min ) then
                dista_min = dista
                ipoin_min = ipoin
             end if
          end if
       end do

       inout_test(kpoin) = ipoin_min
       dimin_test(kpoin) = dista_min

    end do

  end subroutine COU_WET_POINTS_NEAREST_BOUNDARY_NODE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          check which of them host them
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_NEAREST_ELEMENT_NODE(npoin_test,coord_test,dimin_test,inout_test,lnsou)
    integer(ip), intent(in)          :: npoin_test                    !< Number of points to test
    real(rp),    intent(in)          :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    real(rp),    intent(out)         :: dimin_test(npoin_test)        !< Minimum distance to a boundary node
    integer(ip), intent(out)         :: inout_test(npoin_test)        !< Node with minimum distance
    logical(lg), intent(in), pointer :: lnsou(:)
    integer(ip)                      :: ipoin,ibono,ipoin_min
    integer(ip)                      :: idime,kpoin
    real(rp)                         :: dista,dista_min

    do kpoin = 1,npoin_test

       dista_min = huge(1.0_rp)
       ipoin_min = 0

       do ipoin = 1,npoin
          if( lnsou(ipoin) ) then
             dista = dot_product(coord_test(1:ndime,kpoin)-coord(1:ndime,ipoin),coord_test(1:ndime,kpoin)-coord(1:ndime,ipoin))
             if( dista < dista_min ) then
                dista_min = dista
                ipoin_min = ipoin
             end if
          end if
       end do

       inout_test(kpoin) = ipoin_min
       dimin_test(kpoin) = dista_min

    end do

  end subroutine COU_WET_POINTS_NEAREST_ELEMENT_NODE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          Look for a node with the same coordinates
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_SAME_COORDINATE(npoin_test,coord_test,inout_test,lnsou)
    integer(ip),          intent(in)   :: npoin_test                    !< Number of points to test
    real(rp),             intent(in)   :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    integer(ip),          intent(out)  :: inout_test(npoin_test)        !< Node with minimum distance
    logical(lg), pointer, intent(in)   :: lnsou(:)
    integer(ip)                        :: ipoin,ibono
    integer(ip)                        :: idime,kpoin

    if( ndime == 2 ) then

       if( associated(lnsou) ) then
          do kpoin = 1,npoin_test            
             inout_test(kpoin) = 0            
             loop_ipoin_2d_1: do ipoin = 1,npoin
                if( lnsou(ipoin) ) then
                   if( abs(coord_test(1,kpoin)-coord(1,ipoin)) < zeror ) then
                      if( abs(coord_test(2,kpoin)-coord(2,ipoin)) < zeror ) then
                         
                         !if( abs(coord_test(1,kpoin)-0.25000000000000000_rp)<1.0e-6.and.       &
                         !     abs(coord_test(2,kpoin)-0.80000000000000004_rp)<1.0e-6) then
                         !   print*,'a=',coord(1,ipoin),coord(2,ipoin)
                         !   call runend('O.K.!')
                         !end if

                         inout_test(kpoin) = ipoin
                         exit loop_ipoin_2d_1
                      end if
                   end if
                end if
             end do loop_ipoin_2d_1
          end do
       else
          do kpoin = 1,npoin_test            
             inout_test(kpoin) = 0            
             loop_ipoin_2d_2: do ipoin = 1,npoin
                if( abs(coord_test(1,kpoin)-coord(1,ipoin)) < zeror ) then
                   if( abs(coord_test(2,kpoin)-coord(2,ipoin)) < zeror ) then
                      inout_test(kpoin) = ipoin
                      exit loop_ipoin_2d_2
                   end if
                end if
             end do loop_ipoin_2d_2    
          end do
       end if
       
    else
       
       if( associated(lnsou) ) then
          do kpoin = 1,npoin_test          
             inout_test(kpoin) = 0
             loop_ipoin_3d_1: do ipoin = 1,npoin
                if( lnsou(ipoin) ) then
                   if( abs(coord_test(1,kpoin)-coord(1,ipoin)) < zeror ) then
                      if( abs(coord_test(2,kpoin)-coord(2,ipoin)) < zeror ) then
                         if( abs(coord_test(3,kpoin)-coord(3,ipoin)) < zeror ) then
                            inout_test(kpoin) = ipoin
                            exit loop_ipoin_3d_1
                         end if
                      end if
                   end if
                end if
             end do loop_ipoin_3d_1
          end do
       else
          do kpoin = 1,npoin_test          
             inout_test(kpoin) = 0
             loop_ipoin_3d_2: do ipoin = 1,npoin
                if( abs(coord_test(1,kpoin)-coord(1,ipoin)) < zeror ) then
                   if( abs(coord_test(2,kpoin)-coord(2,ipoin)) < zeror ) then
                      if( abs(coord_test(3,kpoin)-coord(3,ipoin)) < zeror ) then
                         inout_test(kpoin) = ipoin
                         exit loop_ipoin_3d_2
                      end if
                   end if
                end if
             end do loop_ipoin_3d_2
          end do
       end if
       
    end if
    
  end subroutine COU_WET_POINTS_SAME_COORDINATE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          check which of them host them
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_HOST_BOUNDARIES(&
       npoin_test,coord_test,dimin_test,inout_test,shapf_test,kdtree_source,lbsou)
    integer(ip),                        intent(in)    :: npoin_test                    !< Number of points to test
    real(rp),                           intent(in)    :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    real(rp),                           intent(out)   :: dimin_test(npoin_test)        !< Minimum distance to a boundary node
    integer(ip),                        intent(out)   :: inout_test(npoin_test)        !< Node with minimum distance
    real(rp),                           intent(out)   :: shapf_test(mnodb,npoin_test)  !< Shape function in elements
    type(typ_kdtree),                   intent(inout) :: kdtree_source                 !< KD-tree
    logical(lg),      pointer, optional,intent(in)    :: lbsou(:)                      !< List of source boundaries
    integer(ip)                                       :: ipoin,iboun_min
    integer(ip)                                       :: kpoin,inodb
    integer(ip)                                       :: pnodb,pblty,ifoun,ptopo
    real(rp)                                          :: dista_min,proje(3)
    real(rp)                                          :: deriv(3*mnode)
    real(rp)                                          :: coloc(3),toler
    real(rp)                                          :: bocod(ndime,mnodb)
    logical(lg)                                       :: if_mask

    if_mask = .false.
    if( present(lbsou) ) then
       if( associated(lbsou) ) if_mask = .true.
    end if
    
    if( kdtree_source % nboun /= 0 ) then

       toler = 0.01_rp

       do kpoin = 1,npoin_test 

          dista_min = huge(1.0_rp)
          !
          ! KDTREE returns the boundary number in global numbering because KDTree was constructed 
          ! using a permutation array
          !
          if( if_mask ) then
             call kdtree_nearest_boundary(coord_test(1,kpoin),kdtree_source,iboun_min,dista_min,proje,LMASK=lbsou)
          else
             call kdtree_nearest_boundary(coord_test(1,kpoin),kdtree_source,iboun_min,dista_min,proje)
          end if
          ! if( dista_min < 0 ) print*, "DEBUG: dista min= ",dista_min,coord_test(:,kpoin),kpoin
          dista_min = abs(dista_min)          
          ifoun     = 0
          pblty     = abs(ltypb_cou(iboun_min))
          pnodb     = lnnob_cou(iboun_min)   
          ptopo     = ltopo(pblty)
          do inodb = 1,pnodb
             ipoin = lnodb_cou(inodb,iboun_min)
             bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
          end do
          !print*,ndime,pblty,ptopo,pnodb

          call elmgeo_natural_coordinates_on_boundaries(&
               ndime,pblty,pnodb,bocod,         &
               shapf_test(1,kpoin),deriv,proje, & 
               coloc,ifoun,toler)
          !
          ! To avoid problems temporarly
          !

          ifoun = 1_ip

          !
          !if(ifoun/=0) write(*,'(i2,7(1x,e12.6))') kfl_paral,coord_test(1:2,kpoin),kdtree_typ % coord(1:ndime,kdtree_typ % lnodb_cou1,iboun_min)),kdtree_typ % coord(1:ndime,kdtree_typ % lnodb_cou2,iboun_min))
          !if(ifoun/=0) write(*,'(7(1x,e12.6))') kdtree_typ % coord(1:ndime,kdtree_typ % lnodb_cou1,iboun_min)),coord(1:ndime,lnodb_cou1,iboun_min))

          if( ifoun == 0 ) then
             !write(6,'(a,3(1x,i5))')             'A=',kfl_paral,kpoin,iboun_min,nboun
             !write(6,'(a,2(1x,i4),3(1x,e12.6))') 'B=',kfl_paral,kpoin,coloc

             !write(6,'(a,i2,7(1x,i5))')          'C=',kfl_paral,kpoin,lnodb_cou1,iboun_min),lnodb_cou2,iboun_min),lnodb_cou3,iboun_min),lnodb_cou4,iboun_min)
             !write(6,'(a,i2,i5,7(1x,e12.6))')    'D=',kfl_paral,kpoin,coord_test(1:3,kpoin)
             !write(6,*) '1: ',coord(1:3,lnodb_cou1,iboun_min))
             !write(6,*) '2: ',coord(1:3,lnodb_cou2,iboun_min))
             !write(6,*) '3: ',coord(1:3,lnodb_cou3,iboun_min))
             !write(6,*) '4: ',coord(1:3,lnodb_cou4,iboun_min))
             !!coord(1:ndime,kdtree_typ % lnodb_cou1,iboun_min)),kdtree_typ % coord(1:ndime,kdtree_typ % lnodb_cou2,iboun_min))
             !call iofile_flush_unit(6_ip)
             !call runend('O.K.!')
             inout_test(kpoin) = 0
             dimin_test(kpoin) = dista_min
          else
             inout_test(kpoin) = iboun_min
             dimin_test(kpoin) = dista_min    
          end if

       end do

    end if

  end subroutine COU_WET_POINTS_HOST_BOUNDARIES

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    11/03/2014
  !> @brief   Check in which elements are the test points
  !> @details For each partition whcih have received some points,
  !>          check which of them host them
  !
  !----------------------------------------------------------------------

  subroutine COU_WET_POINTS_HOST_BOUNDARY_VECTOR(&
       npoin_test,coord_test,dimin_test,inout_test,shapf_test,kdtree_source)
    integer(ip),                        intent(in)    :: npoin_test                    !< Number of points to test
    real(rp),                           intent(in)    :: coord_test(ndime,npoin_test)  !< Coordinates of test points
    real(rp),                           intent(out)   :: dimin_test(npoin_test)        !< Minimum distance to a boundary node
    integer(ip),                        intent(out)   :: inout_test(npoin_test)        !< Node with minimum distance
    real(rp),                           intent(out)   :: shapf_test(mnodb,npoin_test)  !< Shape function in elements
    type(typ_kdtree),                   intent(inout) :: kdtree_source                 !< KD-tree
    integer(ip)                                       :: ipoin,iboun_min
    integer(ip)                                       :: kpoin,inodb
    integer(ip)                                       :: pnodb,pblty,ifoun,ptopo
    real(rp)                                          :: dista_min,proje(3)
    real(rp)                                          :: deriv(3*mnode)
    real(rp)                                          :: coloc(3),toler
    real(rp)                                          :: bocod(ndime,mnodb)


    call runend('COU_WET_POINTS_HOST_BOUNDARY_VECTOR is not ready')

    if( kdtree_source % nboun /= 0 ) then

       toler = 0.01_rp

       do kpoin = 1,npoin_test 
          !
          ! Look for boundary for wet point coord_test(:,kpoin)
          !
          ! shapf_test(1:mnodb,kpoin)
          ! inout_test(kpoin) = iboun_min
          ! dimin_test(kpoin) = dista_min
          !call elmgeo_natural_coordinates(      &
          !     ndime,pblty,ptopo,pnodb,bocod,   &
          !     shapf_test(1,kpoin),deriv,proje, &
          !     coloc,ifoun,toler)
       end do

    end if

  end subroutine COU_WET_POINTS_HOST_BOUNDARY_VECTOR

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    23/09/2014
  !> @brief   List of wet nodes
  !> @details Give the list of source node for a specific coupling
  !>          of for all couplings
  !
  !----------------------------------------------------------------------

  subroutine COU_LIST_SOURCE_NODES(list_source_nodes,kcoup,what_to_do)
    logical(lg),  intent(inout), pointer  :: list_source_nodes(:)
    integer(ip),  intent(in),    optional :: kcoup
    character(*), intent(in),    optional :: what_to_do
    integer(ip)                           :: icoup_ini,icoup_end,icoup
    integer(ip)                           :: kpoin,ipoin
    !
    ! Allocate of necessary
    !
    if( associated(list_source_nodes) ) then
       if( present(what_to_do) ) then
          if( trim(what_to_do) == 'INITIALIZE' ) then
             do ipoin = 1,size(list_source_nodes)
                list_source_nodes(ipoin) = .false.
             end do
          end if
       end if
    else
       call memory_alloca(memor_cou,'COU_LIST_SOURCE_NODES','list_source_nodes',list_source_nodes,npoin)
    end if
    !
    ! Bounds
    !
    if( present(kcoup) ) then
       if( kcoup /= 0 ) then
          icoup_ini = kcoup
          icoup_end = kcoup
       else
          icoup_ini = 1
          icoup_end = mcoup          
       end if
    else
       icoup_ini = 1
       icoup_end = mcoup
    end if
    !
    ! Activate node
    !
    do icoup = icoup_ini,icoup_end
       do kpoin = 1,coupling_type(icoup) % geome % npoin_source
          ipoin = coupling_type(icoup)  % geome % lpoin_source(kpoin)
          list_source_nodes(ipoin) = .true.
       end do
    end do

  end subroutine COU_LIST_SOURCE_NODES

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   If there is at least one zone coupling
  !> @details If there is at least one zone coupling
  !
  !----------------------------------------------------------------------

  function THERE_EXISTS_A_ZONE_COUPLING()
    integer(ip) :: icoup
    logical(lg) :: THERE_EXISTS_A_ZONE_COUPLING

    THERE_EXISTS_A_ZONE_COUPLING = .false.
    do icoup = 1,mcoup
       if( coupling_type(icoup) % zone_target + coupling_type(icoup) % zone_source /= 0 ) then
          THERE_EXISTS_A_ZONE_COUPLING = .true.
          return
       end if
    end do
  end function THERE_EXISTS_A_ZONE_COUPLING

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/03/2014
  !> @brief   Initialize a value on target nodes
  !> @details Initialize a value on target nodes only if the coupling
  !>          is not on whole mesh
  !>
  !----------------------------------------------------------------------

  subroutine COU_PUT_VALUE_ON_TARGET_IP_1(value_in,xarray,kcoup)
    integer(ip), intent(in)              :: value_in
    integer(ip), intent(inout), pointer  :: xarray(:)
    integer(ip), intent(in),    optional :: kcoup
    integer(ip)                          :: icoup_ini,icoup_fin
    integer(ip)                          :: kpoin,ipoin,icoup

    if( mcoup > 0 ) then 
       if( present(kcoup) ) then
          icoup_ini = kcoup
          icoup_fin = kcoup
       else
          icoup_ini = 1
          icoup_fin = mcoup
       end if

       do icoup = icoup_ini,icoup_fin
          if( coupling_type(icoup) % where_type /= ON_WHOLE_MESH ) then
             do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                if( ipoin > 0 ) xarray(ipoin) = value_in
             end do
          end if
       end do
    end if

  end subroutine COU_PUT_VALUE_ON_TARGET_IP_1

  subroutine COU_PUT_VALUE_ON_TARGET_IP_2(value_in,xarray,kcoup)
    integer(ip), intent(in)              :: value_in(*)
    integer(ip), intent(inout), pointer  :: xarray(:,:)
    integer(ip), intent(in),    optional :: kcoup
    integer(ip)                          :: icoup_ini,icoup_fin
    integer(ip)                          :: ndofn,kpoin,ipoin,icoup

    if( mcoup > 0 ) then
       ndofn = size(xarray,1)
       if( present(kcoup) ) then
          icoup_ini = kcoup
          icoup_fin = kcoup
       else
          icoup_ini = 1
          icoup_fin = mcoup
       end if

       do icoup = icoup_ini,icoup_fin
          if( coupling_type(icoup) % where_type /= ON_WHOLE_MESH ) then
             do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                if( ipoin > 0 ) xarray(1:ndofn,ipoin) = value_in(1:ndofn)
             end do
          end if
       end do
    end if

  end subroutine COU_PUT_VALUE_ON_TARGET_IP_2

  subroutine COU_PUT_VALUE_ON_TARGET_IP_12(value_in,xarray,kcoup)
    integer(ip), intent(in)              :: value_in
    integer(ip), intent(inout), pointer  :: xarray(:,:)
    integer(ip), intent(in),    optional :: kcoup
    integer(ip)                          :: icoup_ini,icoup_fin
    integer(ip)                          :: ndofn,kpoin,ipoin,icoup

    if( mcoup > 0 ) then
       ndofn = size(xarray,1)
       if( present(kcoup) ) then
          icoup_ini = kcoup
          icoup_fin = kcoup
       else
          icoup_ini = 1
          icoup_fin = mcoup
       end if

       do icoup = icoup_ini,icoup_fin
          if( coupling_type(icoup) % where_type /= ON_WHOLE_MESH ) then
             do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                if( ipoin > 0 ) xarray(1:ndofn,ipoin) = value_in
             end do
          end if
       end do
    end if

  end subroutine COU_PUT_VALUE_ON_TARGET_IP_12

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    14/10/2014
  !> @brief   Change fixity array
  !> @details Modify fixity array on target:
  !>          - Coupling between zones:
  !>            UNKNOWN type: set fixity to FIXED_UNKNOWN whenener it 
  !>            is different from 0. Can be forced by using "FORCE"
  !>          - Coupling between subdomains:
  !>            RESIDUAL type: free the nodes if FORCE SOBDOMAIN 
  !             is present
  !>
  !----------------------------------------------------------------------

  subroutine COU_SET_FIXITY_ON_TARGET(variable,imodu,kfl_fixno,what)
    character(*), intent(in)              :: variable
    integer(ip),  intent(in)              :: imodu
    integer(ip),  intent(inout), pointer  :: kfl_fixno(:,:)
    character(*), intent(in), optional    :: what
    integer(ip)                           :: ipoin,kpoin,npoin_wet
    integer(ip)                           :: ndofn,kdofn,kpoin_fixed
    integer(ip)                           :: icoup,iffix_max,idofn
    integer(ip)                           :: code_target,zone_target
    logical(lg)                           :: force_subdomain

    if( mcoup == 0 ) return
    !
    ! What to do
    !
    iffix_max       = 0
    ndofn           = size(kfl_fixno,1)
    force_subdomain = .false.
    if( present(what) ) then
       if( trim(what) == 'FORCE ZONE' ) then
          iffix_max = huge(1_ip)
       else if( trim(what) == 'FREE BETWEEN SUBDOMAINS' ) then
          force_subdomain = .true.
       end  if
    end if

    do icoup = 1,mcoup

       if( coupling_type(icoup) % where_type /= ON_WHOLE_MESH ) then

          code_target  = coupling_type(icoup) % code_target
          zone_target  = coupling_type(icoup) % zone_target
          color_target = coupling_type(icoup) % color_target
          npoin_wet    = coupling_type(icoup) % wet % npoin_wet

          if( I_AM_IN_COLOR(color_target) ) then

             if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
                !
                ! Between subdomains: just check not all dofs are prescribed
                !
                if( force_subdomain ) then
                   kpoin_fixed = 0
                   do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                      ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                      kfl_fixno(1:ndofn,ipoin) = 0
                   end do
                else
                   kpoin_fixed = 0
                   do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                      ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                      if( minval(kfl_fixno(1:ndofn,ipoin)) > 0 ) kpoin_fixed = kpoin_fixed + 1
                   end do
                   if( kpoin_fixed == npoin_wet ) then
                      kpoin = 1
                   else
                      kpoin = 0
                   end if
                   call PAR_MIN(kpoin,'IN CURRENT TARGET COLOR')
                   if( kpoin == 1 ) then
                      call runend('ALL DOFS ARE PRESCIBED ON INTERFACE')
                   end if
                end if

             else if( coupling_type(icoup) % kind == BETWEEN_ZONES .and. lzone(imodu) == zone_target  ) then
                !
                ! Between zones
                !
                ! Unknown type:  force iffix
                ! Residual type: just check not all dofs are prescribed
                !
                if( coupling_type(icoup) % variable == variable(1:5) ) then

                   if( coupling_type(icoup) % what == UNKNOWN  ) then
                      do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                         ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                         do idofn = 1,ndofn
                            if( kfl_fixno(idofn,ipoin) <= iffix_max ) kfl_fixno(idofn,ipoin) = FIXED_UNKNOWN
                         end do
                      end do

                   else if( coupling_type(icoup) % what == RESIDUAL ) then

                      kpoin_fixed = 0
                      do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                         ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                         kdofn = 0
                         if( minval(kfl_fixno(1:ndofn,ipoin)) > 0 ) kpoin_fixed = kpoin_fixed + 1
                      end do
                      if( kpoin_fixed == npoin_wet ) then
                         kpoin = 1
                      else
                         kpoin = 0
                      end if
                      call PAR_MIN(kpoin,'IN CURRENT TARGET COLOR')
                      if( kpoin == 1 ) then
                         print*,'popopo=',current_code,' ',variable
                         call runend('ALL DOFS ARE PRESCIBED ON INTERFACE')
                      end if

                   end if

                end if

             end if

          end if

       end if

    end do

  end subroutine COU_SET_FIXITY_ON_TARGET

  !----------------------------------------------------------------------
  !>
  !> @author  J.C. Cajas
  !> @date    02/06/2016
  !> @brief   Predict interface values
  !> @details Predict coupling values when advancing in time
  !>          
  !>
  !----------------------------------------------------------------------
  subroutine COU_TEMPORAL_PREDICTOR(icoup)

    integer(ip), intent(in) :: icoup
    integer(ip)             :: idime, kpoin, npoin_wet


    if( coupling_type(icoup) % temporal_predictor /= 1_ip ) return

    npoin_wet = coupling_type(icoup) % wet % npoin_wet

    if( .not. associated( coupling_type(icoup) % values_converged ) .and. npoin_wet > 0 ) &
         &call runend(' Sure you made two time steps before making a prediction? see cou_update_points_values, runend from cou_temporal_predictor ')
    !
    ! Call madame Sazu and make your prediction, it will be stored in 
    ! coupling_type(icoup) % values_converged(:,:,1) in order to reuse the array
    !
    ! prediction(k+1) = 5/2 * values_conv(k) - 2 * values_conv(k-1) + 1/2 * values_conv(k-2)
    !

    if( coupling_type(icoup) % temporal_predictor_order == 0_ip ) then
      return ! Value is the one of the last time step. Saved in convergence subroutine

    elseif( coupling_type(icoup) % temporal_predictor_order == 1_ip ) then
      call runend('COU_TEMPORAL_PREDICTOR: first order predictor not coded')
    
    elseif( coupling_type(icoup) % temporal_predictor_order == 2_ip ) then
      do kpoin = 1_ip, npoin_wet
        do idime = 1_ip, ndime
          coupling_type(icoup) % values_converged(idime,kpoin,1_ip) = 2.5_rp * coupling_type(icoup) % values_converged(idime,kpoin,1_ip) - &
               & 2.0_rp * coupling_type(icoup) % values_converged(idime,kpoin,2_ip) + 0.5_rp * coupling_type(icoup) % values_converged(idime,kpoin,3_ip)
        end do
      end do

    else
      call runend('MOD_COUPLINGS: temporal predictor order not programmed')

    endif

  end subroutine COU_TEMPORAL_PREDICTOR


      !----------------------------------------------------------------------
      !>
      !> @author  J.C. Cajas
      !> @date    02/06/2016
      !> @brief   Q-N Broyden ('bad') approximation 
      !> @details Broyden (bad) Q-N approximation
      !>          
      !>
      !----------------------------------------------------------------------
  subroutine COU_BROYDEN_BAD(xxnew_o,coupli_o)

    real(rp),      pointer,   intent(inout) :: xxnew_o(:,:)
    type(typ_color_coupling), intent(inout) :: coupli_o
    integer(ip)                             :: ndofn,ipoin,npoin_wet,npoin_wet_total
    integer(ip)                             :: idofn,kpoin,jpoin
    integer(ip)                             :: initial_index, itime
    integer(ip)                             :: ntime_wet,ipart
    integer(ip)                             :: target_rank, target_size

    ! Broyden auxiliary variables
    real(rp), pointer                       :: deltaf(:), deltag(:), deltav(:), auxjac(:,:)
    real(rp), pointer                       :: my_dex(:), my_def(:), my_dev(:)
    real(rp)                                :: normv, relax

    integer(ip), pointer                    :: npoin_wet_all(:)


    nullify(my_dex)
    nullify(my_def)
    nullify(my_dev)

    nullify(deltaf)
    nullify(deltav)

    nullify(deltag)
    nullify(auxjac)

    nullify(npoin_wet_all)

    if( INOTMASTER ) then
       !
       ! Degrees of freedom
       !
       ndofn     = size(xxnew_o,1)
       npoin_wet = coupli_o % wet % npoin_wet
       ntime_wet = 3

    else

       ndofn     = 0
       npoin_wet = 0
       ntime_wet = 0

    end if
    !
    ! Size of the whole wet surface
    !    
    if( .not. associated(npoin_wet_all) )then

       call PAR_COMM_RANK_AND_SIZE(target_rank,target_size,'IN CURRENT TARGET COLOR')
       allocate(npoin_wet_all(0:target_size-1))
       !call memory_alloca(memor_cou,'MOD_COUPLINGS',' npoin_all ',npoin_wet_all,target_size)

       call PAR_ALLGATHER(npoin_wet,npoin_wet_all,1_4,'IN CURRENT TARGET COLOR')

       npoin_wet_total = 0_rp
       do ipart = 0_ip, target_size-1_ip
          if( ipart == target_rank )initial_index = npoin_wet_total
          npoin_wet_total = npoin_wet_total + npoin_wet_all(ipart)
       end do

       npoin_wet_all = ndime * npoin_wet_all


    end if
    ! print*, "DEBUG: npoin_wet ", npoin_wet_all, target_rank, npoin_wet_total
    !
    ! Memory allocation
    !
    if( INOTMASTER )then

       if( .not. associated(coupli_o % values) )          &
            call memory_alloca(memor_cou,'MOD_COUPLINGS','values',coupli_o % values,ndofn,npoin_wet,3_ip)
       if( .not. associated(coupli_o % values_predicted) )&
            call memory_alloca(memor_cou,'MOD_COUPLINGS','values_predicted',coupli_o % values_predicted,ndofn,npoin_wet)

       if( .not. associated(auxjac) )                     &
            call memory_alloca(memor_cou,'MOD_COUPLINGS',' auxjac ', auxjac ,ndofn * npoin_wet, ndofn * npoin_wet_total)
       if( .not. associated(coupli_o % jacobian_inverse) )&
            call memory_alloca(memor_cou,'MOD_COUPLINGS','jacobian_inverse',coupli_o % jacobian_inverse,ndofn,npoin_wet,npoin_wet_total)

       if( .not. associated(my_dex) )                     &
            call memory_alloca(memor_cou,'MOD_COUPLINGS',' my_dex ', my_dex ,ndofn * npoin_wet)
       if( .not. associated(my_def) )                     &
            call memory_alloca(memor_cou,'MOD_COUPLINGS',' my_def ', my_def ,ndofn * npoin_wet)
       if( .not. associated(deltaf) )                     &
            call memory_alloca(memor_cou,'MOD_COUPLINGS',' deltaf ', deltaf ,ndofn * npoin_wet_total)

       if( .not. associated(my_dev) )                     &
            call memory_alloca(memor_cou,'MOD_COUPLINGS',' my_dev ', my_dev ,ndofn * npoin_wet)
       if( .not. associated(deltav) )                     &
            call memory_alloca(memor_cou,'MOD_COUPLINGS',' deltav ', deltav ,ndofn * npoin_wet_total)

       if( .not. associated(deltag) )                     &
            call memory_alloca(memor_cou,'MOD_COUPLINGS',' deltag ', deltag ,ndofn * npoin_wet)

    else

       if( .not. associated(coupli_o % values) )             call memory_alloca_min( coupli_o % values )
       if( .not. associated(coupli_o % values_predicted) )   call memory_alloca_min( coupli_o % values_predicted )
       if( .not. associated(coupli_o % jacobian_inverse) )   call memory_alloca_min( coupli_o % jacobian_inverse )


       call memory_alloca_min( my_dex )
       call memory_alloca_min( my_def )
       call memory_alloca_min( my_dev )

       if( .not. associated(deltaf) )                     &
            call memory_alloca(memor_cou,'MOD_COUPLINGS',' deltaf ', deltaf ,ndime * npoin_wet_total)
       if( .not. associated(deltav) )                     &
            call memory_alloca(memor_cou,'MOD_COUPLINGS',' deltav ', deltav ,ndime * npoin_wet_total)


       call memory_alloca_min( deltag )
       call memory_alloca_min( auxjac )

       npoin_wet_total = 0_rp

    end if
    !
    ! Save old relaxed values
    !
    do itime = ntime_wet,2,-1
       do ipoin = 1,npoin_wet
          do idofn = 1,ndofn
             coupli_o % values(idofn,ipoin,itime) = coupli_o % values(idofn,ipoin,itime-1)
          end do
       end do
    end do
    !
    ! Search the roots of f(x)-x
    !
    do ipoin = 1_ip, npoin_wet
       do idofn = 1_ip, ndofn
          xxnew_o(idofn,ipoin)                 = xxnew_o(idofn,ipoin) - coupli_o % values(idofn,ipoin,2)
          my_dev( (ipoin-1_ip)*ndofn + idofn ) = xxnew_o(idofn,ipoin)
       end do
    end do
    !
    ! First iteration is performed with an initial guess of the Jacobian inverse
    !
    if( coupling_driver_iteration(1_ip) < 2_ip ) then

       relax = coupli_o % relax

       if( coupling_driver_iteration(1_ip) == 1_ip ) then

          do ipoin = 1_ip,npoin_wet
             do jpoin = 1_ip, npoin_wet_total
                do idofn = 1_ip,ndofn
                   coupli_o % jacobian_inverse(idofn,ipoin,jpoin) = 0_rp
                end do
             end do
          end do

          do ipoin = 1_ip,npoin_wet
             do idofn = 1_ip,ndofn
                coupli_o % jacobian_inverse(idofn,ipoin,initial_index+ipoin) =-relax
             end do
          end do

       end if

       do ipoin = 1_ip, npoin_wet
          do idofn = 1_ip, ndofn

             my_dex( (ipoin-1_ip)*ndofn + idofn ) = coupli_o % values(idofn,ipoin,2_ip) - coupli_o % values(idofn,ipoin,3_ip)
             my_def( (ipoin-1_ip)*ndofn + idofn ) = xxnew_o(idofn,ipoin) - coupli_o % values_predicted(idofn,ipoin)

          end do
       end do

       call PAR_ALLGATHERV(my_dev,deltav,npoin_wet_all,'IN CURRENT TARGET COLOR')

    else

       ! print*, "DEBUG: mas de dos iter ", PAR_MY_CODE_RANK
       !
       ! Auxiliary vectors to compute the correction of the inverse jacobian
       !
       normv  = 0_rp
       do ipoin = 1_ip, npoin_wet_total
          deltaf(ipoin) = 0_rp
          deltav(ipoin) = 0_rp
       end do
       do ipoin = 1_ip, npoin_wet
          deltag(ipoin) = 0_rp
       end do

       do ipoin = 1_ip, npoin_wet
          do idofn = 1_ip, ndofn
             my_dex( (ipoin-1_ip)*ndofn + idofn ) = coupli_o % values(idofn,ipoin,2_ip) - coupli_o % values(idofn,ipoin,3_ip)
             my_def( (ipoin-1_ip)*ndofn + idofn ) = xxnew_o(idofn,ipoin) - coupli_o % values_predicted(idofn,ipoin)

             ! deltaf( (initial_index+ipoin-1_ip)*ndofn + idofn ) =  my_def( (ipoin-1_ip)*ndofn + idofn )
             ! deltav( (initial_index+ipoin-1_ip)*ndofn + idofn ) =  my_dev( (ipoin-1_ip)*ndofn + idofn )

          end do
       end do

       call PAR_ALLGATHERV(my_def,deltaf,npoin_wet_all,'IN CURRENT TARGET COLOR')
       call PAR_ALLGATHERV(my_dev,deltav,npoin_wet_all,'IN CURRENT TARGET COLOR')

       do ipoin = 1_ip, npoin_wet_total
          do idofn = 1_ip, ndofn
             normv = normv + deltaf( (ipoin-1_ip)*ndofn + idofn ) * deltaf( (ipoin-1_ip)*ndofn + idofn )
          end do
       end do

       normv = 1_rp / ( normv + zeror )

       do ipoin = 1_ip, npoin_wet
          do jpoin = 1_ip, npoin_wet_total
             do idofn = 1_ip, ndofn
                deltag( (ipoin-1_ip)*ndofn + idofn ) = deltag( (ipoin-1_ip)*ndofn + idofn ) + &
                     & coupli_o % jacobian_inverse(idofn,ipoin,jpoin) * deltaf((jpoin-1_ip)*ndofn + idofn )
             end do
          end do
       end do
       do ipoin = 1_ip, npoin_wet
          do idofn = 1_ip, ndofn
             deltag((ipoin-1_ip)*ndofn + idofn ) = ( my_dex( (ipoin-1_ip)*ndofn + idofn ) - &
                  &deltag( (ipoin-1_ip)*ndofn + idofn ) ) * normv
          end do
       end do
       auxjac = 0_rp
       if( INOTMASTER .and. npoin_wet /= 0 )then
          call maths_outer_product(deltag,deltaf,auxjac)
       end if
       !
       ! Update the Jacobian approximation
       !
       do jpoin = 1_ip, npoin_wet_total
          do ipoin = 1_ip, npoin_wet
             do idofn = 1_ip, ndofn

                coupli_o % jacobian_inverse(idofn,ipoin,jpoin) = coupli_o % jacobian_inverse(idofn,ipoin,jpoin) + &
                     &auxjac( (ipoin-1_ip)*ndofn+idofn ,(jpoin-1_ip)*ndofn+idofn  )

             end do
             ! print*, "DEBUG: jacobian_inverse ", ipoin, jpoin, coupli_o % jacobian_inverse(:,ipoin,jpoin)
          end do
       end do

    end if      ! Iterations greater thar 1
    !
    ! Save unmodified results
    !
    do ipoin = 1_ip, npoin_wet
       do idofn = 1_ip, ndofn
          coupli_o % values_predicted(idofn,ipoin) = xxnew_o(idofn,ipoin)
       end do
    end do
    !
    ! Update of the solution
    !
    do ipoin = 1_ip, npoin_wet
       do idofn = 1_ip, ndofn

          coupli_o % values(idofn,ipoin,1_ip) = 0_rp

       end do
    end do
    do ipoin = 1_ip,npoin_wet
       do jpoin = 1_ip,npoin_wet_total
          do idofn = 1_ip,ndofn 
             coupli_o % values(idofn,ipoin,1_ip) = coupli_o % values(idofn,ipoin,1_ip)+& 
                  coupli_o % jacobian_inverse(idofn,ipoin,jpoin) * deltav( (jpoin-1_ip)*ndofn+idofn )
          end do
       end do
    end do

    do ipoin = 1_ip, npoin_wet
       do idofn = 1_ip, ndofn
          coupli_o % values(idofn,ipoin,1_ip) = coupli_o % values(idofn,ipoin,2_ip) - coupli_o % values(idofn,ipoin,1_ip)
          xxnew_o(idofn,ipoin) = coupli_o % values(idofn,ipoin,1_ip) 
       end do

    end do
    ! print*, "DEBUG: salgo Broyden ", PAR_MY_CODE_RANK

    if( associated(my_def) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltaf ',my_def )
    if( associated(my_dev) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltag ',my_dev )

    ! if( associated(deltax) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltax ',deltax )
    if( associated(deltaf) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltaf ',deltaf )
    if( associated(deltav) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltaf ',deltav )

    if( associated(deltag) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' deltag ',deltag )
    if( associated(auxjac) )call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' auxjac ',auxjac )

    if( associated(npoin_wet_all) ) call memory_deallo( memor_cou,' COU_UPDATE_POINTS_VALUES ',' npoin_wet_all ',npoin_wet_all )

  end subroutine COU_BROYDEN_BAD

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-03-29
  !> @brief   Solver initialization
  !> @details Initialize mask for dot product
  !> 
  !-----------------------------------------------------------------------
  
  subroutine couplings_initialize_solver_mask(solve)

    use def_kintyp, only : soltyp
    use def_coupli, only : ncoup_implicit_n
    use def_coupli, only : ncoup_implicit_d
    use def_coupli, only : mask_cou
    type(soltyp), intent(inout) :: solve          !< Solver structure

    if( ncoup_implicit_n + ncoup_implicit_d > 0 ) then
       solve % kfl_mask =  1
       solve % mask     => mask_cou
    end if
    
  end subroutine couplings_initialize_solver_mask

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-03-29
  !> @brief   Solver initialization
  !> @details Initialize mask for dot product
  !> 
  !-----------------------------------------------------------------------
  
  subroutine couplings_initialize_solver_dirichlet_condition(solve)

    use def_kintyp, only : soltyp
    use def_coupli, only : ncoup_implicit_d
    use def_coupli, only : lcoup_implicit_d
    type(soltyp), intent(inout) :: solve          !< Solver structure
    integer(ip)                 :: kcoup,icoup

    if( INOTMASTER .and. mcoup > 0 ) then
        do kcoup = 1,ncoup_implicit_d
           icoup = lcoup_implicit_d(kcoup)           
        !!!   solve % kfl_dirichlet = 2
        end do
     end if
    
   end subroutine couplings_initialize_solver_dirichlet_condition

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-03-29
  !> @brief   Solver initialization
  !> @details Detect mnodes where reaction is required
  !> 
  !-----------------------------------------------------------------------
  
  subroutine couplings_initialize_solver_reaction(solve)

    use def_kintyp, only : ip,soltyp
    use def_master, only : modul,lzone
    use def_master, only : current_zone
    use def_master, only : mem_modul
    use def_domain, only : npoin,r_dom
    use mod_parall, only : color_target
    use mod_parall, only : color_source
   
    type(soltyp), intent(inout), pointer :: solve(:)          !< Solver structure
    integer(ip)                          :: ireaction
    integer(ip)                          :: num_blocks,jzdom
    integer(ip)                          :: ndofn,icoup,jblok
    integer(ip)                          :: iblok,ipoin
    integer(ip)                          :: ndofn_iblok
    integer(ip)                          :: ndofn_jblok
    integer(ip)                          :: ndofn_block

    ireaction =  0
    ndofn      = solve(1) % ndofn
    num_blocks = solve(1) % num_blocks

    if( solve(1) % block_num == 1 ) then
       !
       ! SOLVE(1) % LPOIN_REACTION(1:NPOIN): Mark the nodes where reaction is required
       ! They are the source nodes
       !
       do icoup = 1,mcoup
          if(    coupling_type(icoup) % what == RESIDUAL      .and. &
               & coupling_type(icoup) % kind == BETWEEN_ZONES ) then

             color_target = coupling_type(icoup) % color_target
             color_source = coupling_type(icoup) % color_source

             if( ireaction == 0 ) then
                if( I_AM_IN_COLOR(color_source) ) then
                   solve(1) % kfl_react = max(solve(1) % kfl_react,1_ip)
                   call memory_alloca(mem_modul(1:2,modul),'SOLVE(1) % LPOIN_REACTION','inivar',solve(1) % lpoin_reaction,npoin)
                end if

                if( I_AM_IN_COLOR(color_target) ) then
                   solve(1) % kfl_bvnat = 1
                   call memory_alloca(mem_modul(1:2,modul),'SOLVE(1) % BVNAT','inivar',solve(1) % block_array(1) % bvnat,ndofn,npoin)
                   solve(1) % bvnat => solve(1) % block_array(1) % bvnat
                   do iblok = 2,num_blocks
                      ndofn_block = solve(1) % block_dimensions(iblok)
                      call memory_alloca(mem_modul(1:2,modul),'SOLVE(1) % BVNAT','inivar',solve(1) % block_array(iblok) % bvnat,ndofn_block,npoin)
                   end do
                end if

             end if

             ireaction = ireaction + 1

             if( INOTMASTER ) call COU_LIST_SOURCE_NODES(solve(1) % lpoin_reaction,icoup)
          end if
       end do
    end if

  end subroutine couplings_initialize_solver_reaction

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-03-29
  !> @brief   Solver initialization
  !> @details Initialize solver according to coupling options
  !> 
  !-----------------------------------------------------------------------
  
  subroutine couplings_initialize_solver()

    use def_kintyp, only : ip,soltyp
    use def_master, only : modul,momod,mmodu
    use def_master, only : current_zone
    use def_master, only : I_AM_IN_ZONE
    use def_master, only : lzone
    use def_master, only : kfl_modul
    use def_solver, only : solve_sol

    integer(ip) :: ivari,jvari

    do modul = 1,mmodu
       current_zone = lzone(modul)
       if( kfl_modul(modul) == 1 .and. associated(momod(modul) % solve) .and. I_AM_IN_ZONE(current_zone) ) then
          do ivari = 1,size(momod(modul) % solve)
             solve_sol => momod(modul) % solve(ivari:)
             call couplings_initialize_solver_mask(solve_sol(1))             
             call couplings_initialize_solver_dirichlet_condition(solve_sol(1))
          end do
          jvari = 1
          do while( jvari <= size(momod(modul) % solve) )
             solve_sol => momod(modul) % solve(jvari:)
             call couplings_initialize_solver_reaction(solve_sol)
             jvari = jvari + momod(modul) % solve(jvari) % block_num
          end do 
       end if
    end do

  end subroutine couplings_initialize_solver
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-02
  !> @brief   Impose the Dirichlet condition
  !> @details Impose the Dirichlet condition on implicit couplings
  !> 
  !-----------------------------------------------------------------------

  subroutine couplings_impose_dirichlet(solve,unkno)

    use def_kintyp, only : soltyp
    use def_coupli, only : ncoup_implicit_d
    use def_coupli, only : lcoup_implicit_d 
    use def_coupli, only : mcoup
  
    type(soltyp), intent(in)    :: solve
    real(rp),     intent(inout) :: unkno(*)
    integer(ip)                 :: icoup,kcoup,ndofn

    ndofn = solve % ndofn
    
    if( INOTMASTER .and. mcoup > 0 ) then
        do kcoup = 1,ncoup_implicit_d
           icoup = lcoup_implicit_d(kcoup)           
           if( solve % kfl_iffix /= 0 ) then
              call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,unkno,mask=solve % kfl_fixno)
           else              
              call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,unkno)
           end if
        end do
     end if
 
   end subroutine couplings_impose_dirichlet
   
   !-----------------------------------------------------------------------
   !> 
   !> @author  houzeaux
   !> @date    2018-05-02
   !> @brief   Check the Dirichlet condition
   !> @details Check if a variable satisfies the Dirichlet condition
   !>          on implicit couplings
   !> 
   !-----------------------------------------------------------------------

   subroutine couplings_check_dirichlet(ndofn,unkno,solve)

     use def_kintyp, only : soltyp
     use def_master, only : lninv_loc
     use def_master, only : intost
     use def_domain, only : npoin
     use def_coupli, only : ncoup_implicit_d
     use def_coupli, only : lcoup_implicit_d 
     use def_coupli, only : mcoup

     integer(ip),  intent(in)             :: ndofn
     real(rp),     intent(inout)          :: unkno(*)
     type(soltyp), intent(in),   optional :: solve
     integer(ip)                          :: icoup,kcoup,ipoin,kpoin,idofn,itotn
     real(rp),     pointer                :: unkno_cpy(:)
     real(rp)                             :: eps
     logical(lg)                          :: mask
     
     if( INOTMASTER .and. mcoup > 0 ) then
        allocate(unkno_cpy(npoin*ndofn))
        unkno_cpy(1:ndofn*npoin) = unkno(1:ndofn*npoin)
        do kcoup = 1,ncoup_implicit_d
           icoup = lcoup_implicit_d(kcoup) 
           call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,unkno)
           do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
              ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
              do idofn = 1,ndofn
                 itotn = (ipoin-1)*ndofn+idofn
                 eps   = abs(unkno_cpy(itotn)-unkno(itotn))
                 mask  = .false.
                 if( present(solve) ) then
                    if( solve % kfl_iffix /= 0 ) then
                       if( solve % kfl_fixno(idofn,ipoin) == 1 ) mask = .true.
                    end if
                 end if
                 if( eps > 1.0e-6_rp .and. .not. mask ) then
                    print*,'WE DO NOT SATISFY THE DIRICHLET CONDITION AT NODE '//intost(lninv_loc(ipoin))//', VALUE =',eps
                 end if
              end do
           end do
        end do
        deallocate(unkno_cpy)
     end if

   end subroutine couplings_check_dirichlet

   subroutine couplings_test_transmission_conditions()
     use def_master
     use def_domain
     use def_coupli
     implicit none
     real(rp), pointer :: xx(:)
     integer(ip)       :: ii,idime,icoup,kcoup,kpoin
     
     allocate(xx(max(1_ip,npoin)))
     do ii = 1,npoin
        xx(ii) = 2.0_rp*coord(2,ii)+3.0_rp
     end do
     
     do kcoup = 1,ncoup_implicit_n
        icoup = lcoup_implicit_n(kcoup)
        do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
           ii = coupling_type(icoup) % wet % lpoin_wet(kpoin)           
           xx(ii) = 0.0_rp
        end do
     end do
     
     do kcoup = 1,ncoup_implicit_n
        icoup = lcoup_implicit_n(kcoup)    
        call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,xx)
        do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
           ii = coupling_type(icoup) % wet % lpoin_wet(kpoin)           
           print*,lninv_loc(ii),xx(ii),abs(xx(ii)-(2.0_rp*coord(2,ii)+3.0_rp))
        end do
     end do
     
   end subroutine couplings_test_transmission_conditions

end module mod_couplings
!> @} 
!-----------------------------------------------------------------------
