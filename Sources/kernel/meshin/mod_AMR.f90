!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_AMR.f90
!> @author  houzeaux
!> @date    2020-03-07
!> @brief   Adaptive mesh refinement
!> @details All tools for adaptive mesh refinement
!-----------------------------------------------------------------------

module mod_AMR

  use def_master
  use def_domain
  use mod_memory
  use mod_communications
  use mod_maths
  use mod_parall
  use def_parall
  use def_coupli
  use mod_couplings
  use mod_mesh_type
  use def_inpout
  use mod_redistribution,           only : commd_npoin
  use mod_redistribution,           only : commd_nelem
  use mod_redistribution,           only : commd_nboun
  use mod_output,                   only : output_open_files
  use mod_ecoute,                   only : ecoute
  use def_kermod,                   only : ndivi
  use def_kintyp_mesh,              only : mesh_type_basic
  use mod_elmgeo,                   only : element_type
  use def_elmtyp,                   only : BOFEM
  use mod_meshes,                   only : meshes_submesh
  use mod_meshes,                   only : meshes_glue_two_meshes
  use mod_meshes,                   only : meshes_list_boundary_elements
  use mod_maths,                    only : maths_heap_sort
  use mod_graphs,                   only : graphs_number_to_linked_list
  use mod_graphs,                   only : graphs_elepoi
  use mod_graphs,                   only : graphs_elepoi_deallocate
  use mod_coupling_memory,          only : cou_initialization
  use mod_coupling_memory,          only : cou_deallocate
  use mod_couplings_communications, only : COU_GENERATE_LOCAL_TRANSMISSION_MATRICES
  use mod_interpolation,            only : COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_messages,                 only : messages_live
  use mod_parall_destructor,        only : parall_destructor
  use mod_par_tools,                only : par_tools_ownership_interfaces
  use mod_par_additional_arrays,    only : par_global_variables_arrays
  use mod_domain,                   only : domain_memory_deallocate 
  use mod_par_global_numbering,     only : par_global_numbering_nodes
  use mod_par_global_numbering,     only : par_global_numbering_elements_boundaries
  use mod_kdtree,                   only : kdtree_construct
  use mod_mpio_par_postpr,          only : posmpio_destructor
  use mod_mpio_par_postpr,          only : posmpio_redistribute
  use mod_mesh_type_basic,          only : mesh_type_basic_valid_mesh
  use mod_mesh_type_basic,          only : mesh_type_basic_collapse_nodes
  use mod_strings,                  only : integer_to_string
  use mod_renumbering,              only : renumbering_node_arrays  
  use mod_error_estimate,           only : error_estimate_mesh_size
  use mod_AMR_interpolate,          only : AMR_interpolate
  use mod_alya2gmsh,                only : alya2gmsh_remeshing
  use mod_mesh_type_basic  
  use def_AMR

  implicit none

  private
  
  character(100)           :: vacal='mod_ADR'
  integer(ip)              :: num_passes
  real(rp),  pointer       :: solut(:)
  !
  ! Subroutines
  !
  public :: AMR
  public :: AMR_initialization
  public :: AMR_parall
  public :: AMR_readat
  
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Initialization
  !> @details AMR initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_variable()

    if( kfl_amr_varia == ID_TEMPE ) then
       solut => tempe(:,1)
    end if
    
  end subroutine AMR_variable
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Initialization
  !> @details AMR initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_initialization()

    nullify(solut)
    num_passes           =  0
    
    num_amr              =  0
    kfl_amr              =  0
    kfl_amr_post         =  2
    kfl_amr_freq         =  1e6
    nelem_amr            =  0
    maxit_amr            =  1
    kfl_amr_varia        =  0
    kfl_amr_background   =  AMR_BACKGROUND_BIN 
    limit_amr_background = 10
    boxes_amr_background = 10
    kfl_size_amr         = 0
    min_size_amr         = 0.0_rp
    max_size_amr         = 0.0_rp

  end subroutine AMR_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Initialization
  !> @details AMR initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_parall()

    call PAR_BROADCAST(kfl_amr)
    call PAR_BROADCAST(kfl_amr_post)
    call PAR_BROADCAST(kfl_amr_freq)
    call PAR_BROADCAST(nelem_amr)
    call PAR_BROADCAST(maxit_amr)
    call PAR_BROADCAST(kfl_amr_varia)
    call PAR_BROADCAST(kfl_amr_background)
    call PAR_BROADCAST(limit_amr_background)
    call PAR_BROADCAST(boxes_amr_background)
    call PAR_BROADCAST(kfl_size_amr)
    call PAR_BROADCAST(min_size_amr)
    call PAR_BROADCAST(max_size_amr)

  end subroutine AMR_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Read data
  !> @details Read input data
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_readat()

    kfl_amr = 1
    call ecoute('ker_readat')
    do while(words(1) /= 'ENDAD' )
       if( words(1) == 'POSTP' ) then
          !
          ! Postprocess file name strategy
          !
          if(      words(2) == 'TAG  ' ) then
             kfl_amr_post = 1
          else if( words(2) == 'DIREC' ) then
             kfl_amr_post = 2
          else if( words(2) == 'NONE ' ) then
             kfl_amr_post = 0
          else
             call runend('PAR_REAPRO: UNKNOWN REPARTITION POSTPROCESS STRATEGY')
          end if
          
       else if( words(1) == 'FREQU' ) then
          !
          ! Partition frequency
          !
          kfl_amr_freq = getint('FREQU',10_ip,'#FREQUENCY FOR REPARTITIONING')
          
       else if( words(1) == 'TARGE' ) then
          !
          ! Partition frequency
          !
          nelem_amr = getint('TARGE',0_ip,'#FREQUENCY FOR REPARTITIONING')
          
       else if( words(1) == 'VARIA' ) then
          !
          ! Partition frequency
          !
          if( words(2) == 'TEMPE' ) then
             kfl_amr_varia = ID_TEMPE
          else
             call runend('MOD_AMR: UNKNOWN VARIABLE')
          end if
          
       else if( words(1) == 'ITERA' ) then
          !
          ! Iterations
          !
          maxit_amr = getint('ITERA',1_ip,'#PARAMETER FOR BACKGROUND MESH')
          
       else if( words(1) == 'BACKG' ) then
          !
          ! Background mesh to prescribe adapted mesh size
          !
          if(      words(2) == 'BIN  ' ) then
             kfl_amr_background = AMR_BACKGROUND_BIN
          else if( words(2) == 'OCTRE' ) then
             kfl_amr_background = AMR_BACKGROUND_OCTREE
          else if( words(2) == 'ORIGI' ) then
             kfl_amr_background = AMR_BACKGROUND_ORIGINAL_MESH
          else if( words(2) == 'OCTBI' ) then
             kfl_amr_background = AMR_BACKGROUND_OCTBIN
          else
             call runend('MOD_AMR: UNKNOWN BACKGROUND MESH TYPE')
          end if
          if( exists('AUTOM') ) then
             limit_amr_background = -1
             boxes_amr_background = -1
          else
             if( exists('MAXIM') ) limit_amr_background = getint('MAXIM',10_ip,'#PARAMETER FOR BACKGROUND MESH OCTREE')
             if( exists('BOXES') ) boxes_amr_background = getint('BOXES',10_ip,'#PARAMETER FOR BACKGROUND MESH BIN')
             if( exists('NUMBE') ) boxes_amr_background = getint('NUMBE',10_ip,'#PARAMETER FOR BACKGROUND MESH BIN')
          end if
          
       else if( words(1) == 'SIZES' ) then
          !
          ! Size strategy
          !
          if( words(2) == 'AUTOM' ) then
             kfl_size_amr = 0
          else if( words(2) == 'PRESC' ) then
             kfl_size_amr = 1
          end if
          if( exists('MINIM') ) min_size_amr = getrea('MINIM',0.0_rp,'#MINIMUM MESH SIZE')
          if( exists('MAXIM') ) max_size_amr = getrea('MAXIM',0.0_rp,'#MAXIMUM MESH SIZE')
          
       end if
       call ecoute('ker_readat')
    end do


  end subroutine AMR_readat
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Main driver
  !> @details Perform iterations
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR()

    type(mesh_type_basic)        :: mesh_new      ! Adapted mesh
    type(mesh_type_basic)        :: mesh_bak      ! Cartesian mesh
    type(mesh_type)              :: mesh_complete ! Complete adapted mesh
    integer(ip),         pointer :: rank_nelem(:)
    real(rp),            pointer :: hh_opt(:)
    real(rp),            pointer :: hh_bak(:)
    integer(ip)                  :: ipass

    num_passes = num_passes + 1

    if( num_passes == kfl_amr_freq ) then
       
       nullify(rank_nelem)
       nullify(hh_opt)
       nullify(hh_bak)
       num_amr = num_amr + 1

       !-----------------------------------------------------------------
       !
       ! Copy original mesh
       !
       !-----------------------------------------------------------------

       call mesh_complete % init('COMPLETE_MESH')
       call mesh_new      % init('INITIAL_MESH')
       call mesh_bak      % init('BACKGROUND_MESH')
       call AMR_copy_mesh(mesh_new)                ! MESH_NEW = MESHE(NDIVI)
       
       call messages_live('---')
       call messages_live('ADAPTIVE MESH REFINEMENT '//integer_to_string(num_amr),'START SECTION')

       !-----------------------------------------------------------------
       !
       ! Compute mesh size
       !
       !-----------------------------------------------------------------

       if( nelem_amr == 0 ) nelem_amr = nelem_total
       call AMR_variable()
       call error_estimate_mesh_size(solut,hh_opt,nelem_opt=nelem_amr)
       call AMR_background_mesh(mesh_new,mesh_bak,hh_opt,hh_bak)

       !-----------------------------------------------------------------
       !
       ! Adaptivity iterations
       !
       !-----------------------------------------------------------------
       
       do ipass = 1,maxit_amr
          call messages_live                     ('ITERATION '//integer_to_string(ipass),'START SECTION')
          call AMR_adapt_mesh                    (mesh_new,mesh_bak,hh_bak)
          call AMR_repartitioning                (mesh_new,rank_nelem,ipass)
          call AMR_new_mesh                      (rank_nelem,mesh_new)
          call mesh_type_basic_parallel_interface(mesh_new)
          call mesh_type_basic_global_numbering  (mesh_new)
          call memory_deallo                     (memor_dom,'RANK_NELEM',trim(vacal),rank_nelem)
          call messages_live                     ('ITERATION','END SECTION')
          !mesh_new % name = 'MESH'//integer_to_string(ipass)
          !call mesh_type_basic_output(mesh_new)
       end do

       call memory_deallo(par_memor,'HH_BAK',vacal,hh_bak)
       call memory_deallo(par_memor,'HH_OPT',vacal,hh_opt)
       call mesh_bak % deallo()

       !-----------------------------------------------------------------
       !
       ! Create new mesh
       !
       !-----------------------------------------------------------------

       !
       ! Create complete mesh from basic one 
       !
       call mesh_type_basic_to_complete(mesh_new,mesh_complete)    ! MESH_COMPLETE => MESH_NEW
       !
       ! Create boundary mesh
       !
       call AMR_boundary_mesh(mesh_complete)
       !
       ! Update domain
       !
       call AMR_update_communication_arrays(mesh_complete)
       call AMR_domain                     (mesh_complete)         ! Missing arrays: LNOCH, KFL_CODNO, etc.
       call AMR_update_original_mesh       (mesh_complete)         ! COORD, LNODS, etc. = MESH % COMPLETE
       call AMR_own_nodes()
       call mesh_type_update_last_mesh()
       !call mesh_new % deallo()
       call mesh_complete % deallo()
       ! 
       ! Destroy
       !
       call memory_allocation_mode('DEALLOCATE BEFORE ALLOCATING')    
       call messages_live('DESTRUCT SOLVERS')
       call solver_destructor()
       call messages_live('DESTRUCT DOMAIN')
       call coupli_destructor()
       call domain_destructor()
       call parall_destructor(COMM_MY_CODE_ARRAY=.false.)
       call posmpio_destructor()
       call output_open_files(AMR=.true.)

       call PAR_DEALLOCATE_COMMUNICATION_ARRAY(commd_npoin,par_memor)
       call PAR_DEALLOCATE_COMMUNICATION_ARRAY(commd_nelem,par_memor)
       call PAR_DEALLOCATE_COMMUNICATION_ARRAY(commd_nboun,par_memor)

       call posmpio_redistribute(FORCE_GENERATE_COMMUNICATOR=.true.)

       call par_global_variables_arrays()
       !
       ! Reconstruct domain and restart run
       !
       call Domain()
       !
       ! Interpolate
       !
       call Interp()    
       !
       ! Deallocate couplings
       !
       call cou_deallocate(coupling_AMR_npoin)     
       call cou_deallocate(coupling_AMR_npoin_nearest) 
       call cou_deallocate(coupling_AMR_nelem)      
       call cou_deallocate(coupling_AMR_nboun)
       !
       ! Beginning of the run
       !
       call Solmem()  
       call Begrun()
       
       call memory_allocation_mode('DO NOT DEALLOCATE BEFORE ALLOCATING')
       
       call messages_live('ADAPTIVE MESH REFINEMENT','END SECTION')
       call messages_live('---')


       num_passes = 0
       
    end if

  end subroutine AMR

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-07
  !> @brief   Copy mesh
  !> @details Copy original mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_copy_mesh(mesh_new)

    type(mesh_type_basic), intent(inout) :: mesh_new
    
    call mesh_new % init     ('INITIAL_MESH')
    call mesh_new % copy     (meshe(ndivi))
    call mesh_new % copy_comm(meshe(ndivi))

  end subroutine AMR_copy_mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-07
  !> @brief   Allocate
  !> @details Allocate new mesh structure
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_mesh_allocate(nelem_loc,npoin_loc,geome,MESH_NAME)

    integer(ip),                 intent(in)    :: nelem_loc
    integer(ip),                 intent(in)    :: npoin_loc
    type(mesh_type_basic),       intent(inout) :: geome
    character(len=*),  optional, intent(in)    :: MESH_NAME
    character(20)                              :: my_mesh_name
    integer(ip)                                :: ii

    if( present(MESH_NAME) ) then
       my_mesh_name = trim(MESH_NAME)
    else
       my_mesh_name = 'MESH'
    end if

    call memory_alloca(memor_dom,trim(my_mesh_name)//' % COORD',trim(vacal),geome % coord,ndime,npoin_loc)
    call memory_alloca(memor_dom,trim(my_mesh_name)//' % LNINV',trim(vacal),geome % lninv_loc,npoin_loc)

    call memory_alloca(memor_dom,trim(my_mesh_name)//' % LNODS',trim(vacal),geome % lnods,mnode,nelem_loc)
    call memory_alloca(memor_dom,trim(my_mesh_name)//' % LTYPE',trim(vacal),geome % ltype,nelem_loc)
    call memory_alloca(memor_dom,trim(my_mesh_name)//' % LEINV',trim(vacal),geome % leinv_loc,nelem_loc)

  end subroutine AMR_mesh_allocate

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Mesh size
  !> @details Transfer mesh size onto a Cartesian mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_background_mesh(mesh_new,mesh_bak,hh_opt,hh_bak)

    type(mesh_type_basic),          intent(inout) :: mesh_new
    type(mesh_type_basic),          intent(inout) :: mesh_bak
    real(rp),              pointer, intent(inout) :: hh_opt(:)
    real(rp),              pointer, intent(inout) :: hh_bak(:)
    real(rp)                                      :: comin(3)
    real(rp)                                      :: comax(3)
    real(rp)                                      :: delta(3),xx(3)
    integer(ip)                                   :: idime,boxes(3)
    integer(ip)                                   :: ipoin,ielem,inode
    integer(ip)                                   :: pnode,pelty
    type(typ_color_coupling)                      :: coupling_AMR
    real(rp),              pointer                :: centroid(:,:)
    type(mesh_type_basic)                         :: mesh_hal

    if( IPARALL ) then
       
       call mesh_hal % init     ('HALO_MESH')
       call mesh_hal % copy     (meshe(ndivi),NPOIN_IN=meshe(ndivi) % npoin_2,NELEM_IN=meshe(ndivi) % nelem_2)

       call mesh_hal % copy_comm(meshe(ndivi))       
       nullify(centroid)
       call mesh_bak % init()
    
       if( mesh_hal % nelem > 0 ) then
          !
          ! Background mesh
          !
          if(      kfl_amr_background == AMR_BACKGROUND_BIN ) then
             !
             ! Bin (Cartesian mesh)
             !
             boxes = (/boxes_amr_background,boxes_amr_background,boxes_amr_background/)
             call mesh_type_basic_bin_mesh(mesh_bak,mesh_hal,boxes,OFFSET=1.0e-6_rp,CENTROID=centroid)

          else if( kfl_amr_background == AMR_BACKGROUND_OCTREE ) then
             !
             ! Octree
             !
             call mesh_type_basic_octree_mesh(mesh_bak,mesh_hal,limit_amr_background,OFFSET=1.0e-6_rp,CENTROID=centroid)
             
          else if( kfl_amr_background == AMR_BACKGROUND_OCTBIN ) then
             !
             ! Bin (Cartesian mesh)
             !
             boxes = (/boxes_amr_background,boxes_amr_background,boxes_amr_background/)
             call mesh_type_basic_octbin_mesh(mesh_bak,mesh_hal,boxes,limit_amr_background,OFFSET=1.0e-6_rp,CENTROID=centroid)

          else if( kfl_amr_background == AMR_BACKGROUND_ORIGINAL_MESH ) then
             !
             ! Original mesh
             !
             mesh_bak = mesh_hal
             allocate(centroid(mesh_bak % ndime,mesh_bak % nelem))
             do ielem = 1,mesh_bak % nelem
                pelty = mesh_bak % ltype(ielem)
                pnode = element_type(pelty) % number_nodes       
                xx    = 0.0_rp
                do inode = 1,pnode
                   ipoin = mesh_bak % lnods(inode,ielem)
                   xx(1:mesh_bak % ndime) = xx(1:mesh_bak % ndime) + mesh_bak % coord(1:mesh_bak % ndime,ipoin)
                end do
                centroid(1:mesh_bak % ndime,ielem) = xx(1:mesh_bak % ndime) / real(pnode,rp)
             end do

          end if

       end if
       call memory_alloca(par_memor,'HH_BAK',vacal,hh_bak,mesh_bak % nelem)

       call cou_initialization(coupling_AMR)

       coupling_AMR % itype                  =  ELEMENT_INTERPOLATION
       coupling_AMR % kfl_multi_source       =  0                          ! Enable multiple source nodes
       coupling_AMR % kfl_lost_wet_points    =  1                          ! Do not crash if wet points are lost
       coupling_AMR % number                 =  2020
       coupling_AMR % color_source           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
       coupling_AMR % color_target           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
       coupling_AMR % commd % PAR_COMM_WORLD =  commd % PAR_COMM_WORLD
       coupling_AMR % wet % npoin_wet        =  mesh_bak % nelem
       coupling_AMR % wet % point_type       =  FLOATING_WET_POINT         ! Wet point is not a node
       coupling_AMR % target_entity          =  FLOATING_TARGET_ENTITY     ! Where coupling is eventuall applied
       coupling_AMR % kfl_source_value       =  VALUES_ON_ELEMENTS
       !
       ! Define wet geometry: select boundary nodes only
       !
       call memory_alloca(par_memor,'COUPLING % WET % COORD_WET' ,'par_interface_exchange',coupling_AMR % wet % coord_wet,ndime,coupling_AMR % wet % npoin_wet)

       do ielem = 1,mesh_bak % nelem
          coupling_AMR % wet % coord_wet(1:mesh_bak % ndime,ielem) = centroid(1:mesh_bak % ndime,ielem)
       end do
       !
       ! Coupling
       !
       call COU_INIT_INTERPOLATE_POINTS_VALUES      (coupling_AMR)
       call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling_AMR)
       call COU_GET_INTERPOLATE_POINTS_VALUES       (hh_opt,hh_bak,coupling_AMR)
       call COU_DEALLOCATE                          (coupling_AMR)     

!!$       block
!!$         !
!!$         ! Postprocess background mesh and centroids
!!$         !
!!$         use def_elmtyp
!!$         use mod_postpr
!!$         type(mesh_type_basic) :: mesh_cen
!!$         integer(ip)           :: ielem
!!$         real(rp),   pointer   :: hh_cen(:)
!!$         call mesh_cen % init  ()
!!$         mesh_cen % npoin =  mesh_bak % nelem
!!$         mesh_cen % nelem =  mesh_bak % nelem
!!$         mesh_cen % ndime =  mesh_bak % ndime
!!$         mesh_cen % mnode =  1
!!$         call mesh_cen % alloca()
!!$         if( mesh_cen % npoin > 0 ) then
!!$            mesh_cen % coord = centroid
!!$            do ielem = 1,mesh_cen % nelem
!!$               mesh_cen % lnods(1,ielem) = ielem
!!$               mesh_cen % ltype(ielem)   = POINT
!!$            end do
!!$         end if
!!$         call mesh_bak % append(mesh_cen)
!!$         call mesh_bak % append(mesh_hal)
!!$         nullify(hh_cen)
!!$         if( mesh_bak % nelem > 0 ) then
!!$            allocate(hh_cen(mesh_bak % nelem))
!!$            hh_cen = 0.0_rp
!!$            do ielem = 1,size(hh_bak)
!!$               hh_cen(ielem) = hh_bak(ielem)
!!$            end do
!!$         end if        
!!$         call mesh_type_basic_parallel(mesh_bak,PAR_COMM_MY_CODE,kfl_paral)
!!$         call mesh_bak % output(filename='back-'//integer_to_string(kfl_paral))
!!$         call mesh_type_basic_output(mesh_bak)
!!$         if( associated(mesh_bak % perme) ) deallocate(mesh_bak % perme)
!!$         if( associated(mesh_bak % permn) ) deallocate(mesh_bak % permn)
!!$         call postpr(hh_cen,(/'SIZES','SCALA','NELEM'/),0_ip,1.0_rp,MESH=mesh_bak)
!!$         call mesh_cen % deallo()
!!$         call runend('O.K.!')
!!$       end block
!!$       
       if( associated(centroid) ) deallocate(centroid)
       call mesh_hal % deallo()

    else if( ISEQUEN ) then

       mesh_bak = mesh_new
       do ielem = 1,mesh_bak % nelem
          hh_bak(ielem) = hh_opt(ielem)
       end do

    end if

  end subroutine AMR_background_mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Determine own nodes
  !> @details Determine own nodes and renumber mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_own_nodes()
    
    integer(ip), pointer :: permr(:)
    integer(ip)          :: ipoin,ii
    
    nullify(permr)
    call par_tools_ownership_interfaces(npoin,commd,permr)    
    call renumbering_node_arrays(permr)

    !block
    !  real(rp), pointer :: tmp(:,:)
    !  allocate(tmp(npoin,2))
    !  if( npoin > 0 ) tmp = tempe
    !  do ipoin = 1,npoin
    !     ii = permr(ipoin)
    !     tempe(ipoin,1) = tmp(ii,1)
    ! end do
    !end block
    !
    ! Own nodes
    !
    npoi1     = commd % npoi1
    npoi2     = commd % npoi2
    npoi3     = commd % npoi3
    npoin_own = commd % npoi3
    call memory_deallo(par_memor,'PERMR',trim(vacal),permr)

  end subroutine AMR_own_nodes
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Copy original mesh
  !> @details Copy current mesh to AMR mesh MESH_NEW
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_copy_original_mesh(mesh_new)

    type(mesh_type_basic), intent(inout) :: mesh_new

    call mesh_new % copy      (meshe(ndivi))
    call mesh_new % copy_comm (meshe(ndivi))

  end subroutine AMR_copy_original_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Repartition mesh
  !> @details Move interface elements to create a new partition
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_repartitioning(mesh_new,rank_nelem,ipass)

    use mod_maths,  only :  maths_sfc_1d_to2d3d_tab
    type(mesh_type_basic),    intent(inout) :: mesh_new
    integer(ip),     pointer, intent(inout) :: rank_nelem(:)
    integer(ip),              intent(in)    :: ipass

    integer(ip),     pointer                :: rank_npoin(:)
    integer(ip),     pointer                :: rank_npoin_new(:)
    integer(ip),     pointer                :: sfc_ord(:)
    integer(ip),     pointer                :: sfc_inv(:)
    integer(ip)                             :: ipoin,ielem,inode
    integer(ip)                             :: min_rank,kpass
    integer(ip)                             :: x,y,z, root
    integer(ip)                             :: ipart, sfc_dim
    integer(ip)                             :: iaux
    real(rp)                                :: raux
    integer(2)                              :: typ, typ_aux

    call messages_live('MOVING INTERFACES')

    if( mesh_new % nelem > 0 ) then
       !
       ! Copy communicator 
       !
       nullify(rank_npoin)
       nullify(rank_npoin_new)
       nullify(sfc_ord)
       nullify(sfc_inv)
       call memory_alloca(memor_dom,'RANK_NPOIN',    trim(vacal),rank_npoin    ,mesh_new % npoin)
       call memory_alloca(memor_dom,'RANK_NPOIN_NEW',trim(vacal),rank_npoin_new,mesh_new % npoin)
       call memory_alloca(memor_dom,'RANK_NELEM',    trim(vacal),rank_nelem    ,mesh_new % nelem)
       !
       ! Define SFC orientation and reordering arrays
       !
       iaux    = mod(7_ip*ipass+num_passes*11_ip,24_ip)
       typ     = int(iaux,kind=2)
       typ_aux = typ
       raux    = real(npart,rp)
       raux    = (raux)**(1.0_rp/3.0_rp)
       root    = ceiling(raux)
       sfc_dim = root**3_ip
       call memory_alloca(memor_dom,'SFC_ORD',       trim(vacal),sfc_ord       ,sfc_dim)
       call memory_alloca(memor_dom,'SFC_INV',       trim(vacal),sfc_inv       ,sfc_dim)    
       do ipart = 1, sfc_dim 
          call  maths_sfc_1d_to2d3d_tab(root,ipart,x,y,z,typ)
          typ = typ_aux
          sfc_ord(ipart)          = (x-1)*root*root+(y-1)*root+z
          sfc_inv(sfc_ord(ipart)) = ipart 
       enddo

       !
       ! Evaluate rank_npoin depending on the SFC ordering
       !
       do ipoin = 1,mesh_new % npoin
          rank_npoin(ipoin)     = sfc_ord(kfl_paral)
          !rank_npoin(ipoin)     = kfl_paral
          rank_npoin_new(ipoin) = huge(1_ip)
       end do

       !
       ! Move interface 
       !
       do kpass = 1,2
          call PAR_INTERFACE_NODE_EXCHANGE(rank_npoin,'MIN',COMM = mesh_new % comm)
          do ielem = 1,mesh_new % nelem
             min_rank = huge(1_ip)
             do inode = 1,nnode(mesh_new % ltype(ielem))
                ipoin    = mesh_new % lnods(inode,ielem)
                min_rank = min(min_rank,rank_npoin(ipoin))
             end do
             do inode = 1,nnode(mesh_new % ltype(ielem))
                ipoin                 = mesh_new % lnods(inode,ielem)
                rank_npoin_new(ipoin) = min(rank_npoin_new(ipoin),min_rank)
             end do
             rank_nelem(ielem) = min_rank
          end do
          do ipoin = 1,mesh_new % npoin
             rank_npoin(ipoin) = rank_npoin_new(ipoin)
          end do
       end do
       !
       ! Revert the SFC reordering
       !
       do ielem = 1,mesh_new % nelem
          rank_nelem(ielem) = sfc_inv(rank_nelem(ielem))
       enddo
      
       call memory_deallo(memor_dom,'RANK_NPOIN',trim(vacal),rank_npoin)
       call memory_deallo(memor_dom,'RANK_NPOIN',trim(vacal),rank_npoin_new)
       call memory_deallo(memor_dom,'SFC_ORD'   ,trim(vacal),sfc_ord)
       call memory_deallo(memor_dom,'SFC_INV'   ,trim(vacal),sfc_inv)

    end if

  end subroutine AMR_repartitioning

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-02
  !> @brief   Adapt mesh
  !> @details Main Adaptation subroutine
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_mesh_size(mesh,mesh_size,mesh_bak)

    type(mesh_type_basic),           intent(in)    :: mesh
    real(rp),              pointer,  intent(inout) :: mesh_size(:)
    type(mesh_type_basic), optional, intent(in)    :: mesh_bak
    integer(ip)                                    :: ielem,ipoin,idime,pnode
    integer(ip)                                    :: inode,pelty
    real(rp)                                       :: xx(3),r

    if( present(mesh_bak) ) then

       nullify(mesh_size)
       allocate(mesh_size(mesh_bak % nelem))
       do ielem = 1,mesh_bak % nelem
          pelty = mesh_bak % ltype(ielem)
          pnode = element_type(pelty) % number_nodes       
          xx = 0.0_rp
          do inode = 1,pnode
             ipoin = mesh_bak % lnods(inode,ielem)
             xx(1:ndime) = xx(1:ndime) + mesh_bak % coord(1:ndime,ipoin)
          end do
          xx = xx / real(pnode,rp)
          mesh_size(ielem) = 0.0_rp
          do idime = 1,ndime
             mesh_size(ielem) = mesh_size(ielem) + sqrt((xx(idime) - 0.5_rp)*(xx(idime) - 0.5_rp))
          end do

          r                = sqrt( (xx(1)-0.5_rp)**2 + (xx(2)-0.5_rp)**2 )
          mesh_size(ielem) = 0.1_rp-0.09_rp*exp(-100.0_rp * (r-0.3_rp)**2)
          !mesh_size(ielem) = 0.05_rp * sqrt(mesh_size(ielem))+0.005_rp
          !mesh_size(ielem) = 0.05_rp
       end do
    else

       nullify(mesh_size)
       allocate(mesh_size(mesh % nelem))
       do ielem = 1,mesh % nelem
          pelty = mesh % ltype(ielem)
          pnode = element_type(pelty) % number_nodes       
          xx = 0.0_rp
          do inode = 1,pnode
             ipoin = mesh % lnods(inode,ielem)
             xx(1:ndime) = xx(1:ndime) + mesh % coord(1:ndime,ipoin)
          end do
          xx = xx / real(pnode,rp)
          mesh_size(ielem) = 0.0_rp
          do idime = 1,ndime
             mesh_size(ielem) = mesh_size(ielem) + sqrt((xx(idime) - 0.5_rp)*(xx(idime) - 0.5_rp))
          end do

          r                = sqrt( (xx(1)-0.5_rp)**2 + (xx(2)-0.5_rp)**2 )
          mesh_size(ielem) = 0.1_rp-0.09_rp*exp(-100.0_rp * (r-0.3_rp)**2)
          !mesh_size(ielem) = 0.05_rp * sqrt(mesh_size(ielem))+0.005_rp
          !mesh_size(ielem) = 0.05_rp
       end do
    end if

  end subroutine AMR_mesh_size

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-02
  !> @brief   Adapt mesh
  !> @details Main Adaptation subroutine
  !>
  !>          MESH_CPY: Original mesh
  !>          MESH_EXT: Extracted mesh with mask
  !>          MESH_NOT: Extracted not valid mesh
  !>          MESH_INT: Interface mesh
  !>          MESH_BOU: Boundary of MESH_EXT
  !>          MESH_NEW: New mesh
  !>
  !>          +--------------------------------+
  !>          |                                |
  !>          |                                |
  !>          |            MESH_CPY            |
  !>          |                                |
  !>          |                                |
  !>          +--------------------------------+
  !>                          ||
  !>                          || Extract interface from mesh 
  !>                          \/
  !>          +---------------------+----------+
  !>          |                     |          |
  !>          |                     |          |
  !>          |      MESH_EXT       | MESH_INT |
  !>          |                     |          |
  !>          |                     |          |
  !>          +---------------------+----------+
  !>                     ||
  !>                     || Extract valid mesh to get manifold mesh
  !>                     \/
  !>          +----------+----------+
  !>          |          |          |
  !>          |          |          |
  !>          | MESH_EXT | MESH_NOT |
  !>          |          |          |
  !>          |          |          |
  !>          +----------+----------+
  !>               ||
  !>               || Extract and append boundary mesh
  !>               \/
  !>          +----------+----------+
  !>          |          |          |
  !>          |          |          |
  !>          | MESH_EXT | MESH_BOU |
  !>          |          |          |
  !>          |          |          |
  !>          +----------+----------+
  !>                    ||
  !>                    ||
  !>                    || REMESHING
  !>                    ||
  !>                    \/
  !>               +----------+
  !>               |          |
  !>               |          |
  !>               | MESH_NEW |
  !>               |          |
  !>               |          |
  !>               +----------+
  !>                    ||
  !>                    || Reconstruct global mesh
  !>                    ||
  !>                    \/
  !>          +----------+----------+----------+
  !>          |          |          |          |
  !>          |          |          |          |
  !>          | MESH_NEW | MESH_INT | MESH_NOT |
  !>          |          |          |          |
  !>          |          |          |          |
  !>          +----------+----------+----------+
  !>                          ||
  !>                          || Append to new mesh
  !>                          ||
  !>                          \/
  !>          +--------------------------------+
  !>          |                                |
  !>          |                                |
  !>          |            MESH_NEW            |
  !>          |                                |
  !>          |                                |
  !>          +--------------------------------+
  !>
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_adapt_mesh(mesh_new,mesh_bak,hh_bak)

    use def_elmtyp
    type(mesh_type_basic),                    intent(inout) :: mesh_new
    type(mesh_type_basic), optional,          intent(in)    :: mesh_bak
    real(rp),              optional, pointer, intent(inout) :: hh_bak(:)
    type(mesh_type_basic)                                   :: mesh_cpy
    type(mesh_type_basic)                                   :: mesh_ext
    type(mesh_type_basic)                                   :: mesh_int
    type(mesh_type_basic)                                   :: mesh_bou
    type(mesh_type_basic)                                   :: mesh_not
    integer(ip)                                             :: ielem,ipoin,idime,pnode,pelty
    integer(ip)                                             :: inode,ii,kk
    integer(ip),                    pointer                 :: invpn(:)
    real(rp),                       pointer                 :: mesh_size(:)
    logical(lg),                    pointer                 :: lmask_nelem(:)
    logical(lg),                    pointer                 :: lmask_npoin(:)
    !
    ! Initialize meshes
    !
    nullify(lmask_npoin,lmask_nelem,invpn)
    
    call mesh_cpy % init('COPY')
    call mesh_ext % init('EXTRACTED_MESH')
    call mesh_int % init('MASK_MESH')
    call mesh_bou % init('BOUNDARY_MESH')
    call mesh_not % init('NON_VALID_MESH')
    
    call messages_live('REMESHING')

    if( mesh_new % nelem > 0 ) then
       !
       ! Copy mesh
       !
       call mesh_cpy % copy     (mesh_new)
       call mesh_cpy % copy_comm(mesh_new)
       call mesh_new % deallo   ()
       call mesh_new % init     ()
       !
       ! Mask elementds with at least one interface node
       !
       call memory_alloca(memor_dom,'LMASK_NELEM',trim(vacal),lmask_nelem,mesh_cpy % nelem)
       call memory_alloca(memor_dom,'LMASK_NPOIN',trim(vacal),lmask_npoin,mesh_cpy % npoin)
       do ielem = 1,mesh_cpy % nelem
          lmask_nelem(ielem) = .true.
       end do
       do ipoin = 1,mesh_cpy % npoin
          lmask_npoin(ipoin) = .true.
       end do

       do ii = 1,mesh_cpy % comm % bound_dim
          ipoin = mesh_cpy % comm % bound_perm(ii)
          lmask_npoin(ipoin) = .false.
       end do
       do ielem = 1,mesh_cpy % nelem
          pelty = mesh_cpy % ltype(ielem)
          do inode = 1,element_type(pelty) % number_nodes
             ipoin = mesh_cpy % lnods(inode,ielem)
             if( .not. lmask_npoin(ipoin) ) lmask_nelem(ielem) = .false.
          end do
       end do
       !
       ! Extract mesh: MESH_EXT = MESH_CPY - MESH_INT
       !
       call mesh_ext % extract(mesh_cpy,lmask_nelem,mesh_cmp=mesh_int)
       !
       ! Create a valid mesh: MESH_EXT = MESH_EXT - MESH_NOT
       !
       call mesh_type_basic_valid_mesh(mesh_ext,mesh_not)    
       !
       ! Mesh size just for testing
       !
       if( present(hh_bak) ) then
          mesh_size => hh_bak
       else
          call AMR_mesh_size(mesh_ext,mesh_size,mesh_bak)
       end if     
       !
       ! Extract and append boundary: MESH_EXT = MESH_EXT + MESH_BOU
       !
       call mesh_bou % boundary_mesh(mesh_ext)
       if( present(mesh_bak) ) then
          call mesh_ext % deallo()
          mesh_ext = mesh_bak
       end if
       call mesh_ext % append(mesh_bou)
       !
       ! Mesh size
       !
       if( kfl_size_amr == 0 .and. associated(mesh_size) ) then
          min_size_amr = minval(mesh_size)
          max_size_amr = maxval(mesh_size)
       end if
       !
       ! Remesh
       !
       call alya2gmsh_remeshing(mesh_new,mesh_ext,mesh_size,min_size_amr,max_size_amr)
       if( mesh_new % npoin > 0 ) mesh_new % lninv_loc = 0
       if( mesh_not % npoin > 0 ) mesh_not % lninv_loc = 0
       if( mesh_new % npoin > 0 ) mesh_new % permn     = 0
       if( mesh_not % npoin > 0 ) mesh_not % permn     = 0
       if( .not. present(hh_bak) ) deallocate(mesh_size)
       !
       ! Append extracted meshes
       !
       call mesh_int % append       (mesh_not)
       call mesh_new % append       (mesh_int)
       call mesh_type_basic_collapse_nodes(mesh_new)
       !
       ! Communication arrays
       !
       call mesh_new % copy_comm(mesh_cpy)       
       
       call memory_alloca(memor_dom,'INVPN',trim(vacal),invpn,mesh_cpy % npoin)
       do ipoin = 1,mesh_new % npoin
          ii       = mesh_new % permn(ipoin)
          if( ii /= 0 ) invpn(ii) = ipoin
       end do
       do ii = 1,mesh_new % comm % bound_dim
          ipoin = mesh_new % comm % bound_perm(ii)
          if( ipoin == 0 ) then
             call runend('THIS SITUATION IS IMPOSSIBLE')
          else
             mesh_new % comm % bound_perm(ii) = invpn(ipoin)
          end if
       end do
       call memory_deallo(memor_dom,'INVPN',trim(vacal),invpn)
       
    end if
    !
    ! Global numbering
    !
    call mesh_type_basic_global_numbering(mesh_new)
    !
    ! Deallocate copy mesh
    !
    ii = mesh_new % nelem
    call PAR_SUM(ii)
    call messages_live('   NUMBER OF ELEMENTS= '//integer_to_string(ii))
   
    call mesh_cpy % deallo()
    call mesh_bou % deallo()
    call mesh_not % deallo()
    call mesh_int % deallo()
    call memory_deallo(memor_dom,'LMASK_NELEM',trim(vacal),lmask_nelem)
    call memory_deallo(memor_dom,'LMASK_NPOIN',trim(vacal),lmask_npoin)

  end subroutine AMR_adapt_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Interface nodes
  !> @details Determine interface nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_new_mesh(rank_nelem,mesh_new)

    integer(ip),           pointer, intent(in)    :: rank_nelem(:)
    type(mesh_type_basic),          intent(inout) :: mesh_new
    integer(ip)                                   :: ielem,ipoin,inode
    integer(ip)                                   :: ineig,jpoin,rank,ipart,ii
    integer(ip)                                   :: kelem,mepoi_loc,ielpo
    integer(4)                                    :: PAR_COMM4
    integer(4)                                    :: SIZE4
    integer(4)                                    :: RANK4
    logical(lg),           pointer                :: list_neighbors(:)
    integer(ip),           pointer                :: nelem_send(:)
    integer(ip),           pointer                :: npoin_send(:)
    integer(ip),           pointer                :: npoin_recv(:)
    integer(ip),           pointer                :: nelem_recv(:)
    type(i1p),             pointer                :: inte_npoin(:)
    type(mesh_type_basic), pointer                :: geome_send(:)
    type(mesh_type_basic), pointer                :: geome_recv(:)
    logical(lg),           pointer                :: lmask(:)
    integer(ip),           pointer                :: lelpo_loc(:)
    integer(ip),           pointer                :: pelpo_loc(:)
    integer(ip),           pointer                :: rank_npoin(:)

    call messages_live('MIGRATING MESH')

    nullify(inte_npoin)
    nullify(list_neighbors)
    nullify(nelem_send)
    nullify(nelem_recv)
    nullify(npoin_send)
    nullify(npoin_recv)
    nullify(geome_send)
    nullify(geome_recv)
    nullify(lmask)
    nullify(pelpo_loc)
    nullify(lelpo_loc)
    nullify(rank_npoin)

    PAR_COMM4 = int(mesh_new % comm % PAR_COMM_WORLD,4)
    SIZE4     = mesh_new % comm % SIZE4
    RANK4     = mesh_new % comm % RANK4

    call memory_alloca(memor_dom,'NELEM_SEND'    ,trim(vacal),nelem_send,npart+1,'INITIALIZE',0_ip)
    call memory_alloca(memor_dom,'NPOIN_SEND'    ,trim(vacal),npoin_send,npart+1,'INITIALIZE',0_ip)
    call memory_alloca(memor_dom,'NELEM_RECV'    ,trim(vacal),nelem_recv,npart+1,'INITIALIZE',0_ip)
    call memory_alloca(memor_dom,'NPOIN_RECV'    ,trim(vacal),npoin_recv,npart+1,'INITIALIZE',0_ip)
    call memory_alloca(memor_dom,'LIST_NEIGHBORS',trim(vacal),list_neighbors,npart)

    !--------------------------------------------------------------------
    !
    ! Determine migrating nodes INTE_NPOIN(IPOIN) % L(:) = RANK
    !
    !--------------------------------------------------------------------

    call graphs_elepoi(mesh_new,mepoi_loc,pelpo_loc,lelpo_loc)
    call memory_alloca(memor_dom,'RANK_NPOIN',trim(vacal),rank_npoin,mepoi)
    call memory_alloca(memor_dom,'INTE_NPOIN',trim(vacal),inte_npoin,mesh_new % npoin)

    do ipoin = 1,mesh_new % npoin
       ipart      = 0
       rank_npoin = 0
       do ielpo = pelpo_loc(ipoin),pelpo_loc(ipoin+1)-1
          ielem = lelpo_loc(ielpo)
          rank  = rank_nelem(ielem)
          if( count( rank_npoin == rank ) == 0 ) then
             ipart                = ipart + 1
             rank_npoin(ipart)    = rank
             list_neighbors(rank) = .true.
          end if
       end do
       call memory_alloca(memor_dom,'INTE_NPOIN % L',trim(vacal),inte_npoin(ipoin) % l,ipart)
       inte_npoin(ipoin) % l(1:ipart) = rank_npoin(1:ipart)
    end do
    call memory_deallo(memor_dom,'RANK_NPOIN',trim(vacal),rank_npoin)
    call graphs_elepoi_deallocate(pelpo_loc,lelpo_loc)
    
    !-----------------------------------------------------------------
    !
    ! Send buffer sizes
    !    
    !-----------------------------------------------------------------
    
    do ipoin = 1,mesh_new % npoin
       do ii = 1,memory_size(inte_npoin(ipoin) % l)
          ipart = inte_npoin(ipoin) % l(ii)
          npoin_send(ipart) = npoin_send(ipart) + 1
       end do
    end do
    do ielem = 1,mesh_new % nelem
       ipart = rank_nelem(ielem)
       nelem_send(ipart) = nelem_send(ipart) + 1
    end do    
    call PAR_ALLTOALL(1_ip,1_ip,npoin_send,npoin_recv,PAR_COMM_IN4=PAR_COMM4)
    call PAR_ALLTOALL(1_ip,1_ip,nelem_send,nelem_recv,PAR_COMM_IN4=PAR_COMM4)
    
    !-----------------------------------------------------------------
    !
    ! Extract submesh GEOME_SEND from MESH_NEW to be sent
    !
    !-----------------------------------------------------------------
    
    allocate(geome_send(npart))
    allocate(geome_recv(npart))
    do ipart = 1,npart
       call geome_send(ipart) % init()
       call geome_recv(ipart) % init()
    end do

    geome_send % mnode = mesh_new % mnode
    geome_send % ndime = mesh_new % ndime
    geome_recv % mnode = mesh_new % mnode
    geome_recv % ndime = mesh_new % ndime

    call memory_alloca(memor_dom,'LMASK',trim(vacal),lmask,mesh_new % nelem)
    do ipart = 1,npart              
       if( nelem_send(ipart) > 0 ) then
          where( rank_nelem == ipart ) lmask = .true.
          call geome_send(ipart) % extract(mesh_new,lmask)
          lmask = .false.
       end if
    end do
    call memory_deallo(memor_dom,'LMASK',trim(vacal),lmask)

    !-----------------------------------------------------------------
    !
    ! Send/Recv submeshes GEOME_SEND/GEOME_RECV
    !    
    !-----------------------------------------------------------------

    do ipart = 1,npart
       call geome_recv(ipart) % alloca(NELEM=nelem_recv(ipart),NPOIN=npoin_recv(ipart))
    end do

    do ipart = 1,npart
       geome_recv(ipart) % nelem = nelem_recv(ipart)
       geome_recv(ipart) % npoin = npoin_recv(ipart)
       call mesh_type_basic_send_recv(geome_send(ipart),geome_recv(ipart),ipart,PAR_COMM4)
    end do
    
    !-----------------------------------------------------------------
    !
    ! Create new mesh MESH_NEW from received meshes GEOME_RECV
    !    
    !-----------------------------------------------------------------

    call mesh_new % deallo()
    call mesh_new % init  ('MESH_NEW')
    
    do ipart = 1,npart
       call mesh_new % merge(geome_recv(ipart))
    end do
    
    !--------------------------------------------------------------------
    !
    ! Deallocate memory
    !
    !--------------------------------------------------------------------

    call memory_deallo(memor_dom,'INTE_NPOIN'    ,trim(vacal),inte_npoin)
    call memory_deallo(memor_dom,'LIST_NEIGHBORS',trim(vacal),list_neighbors)
    call memory_deallo(memor_dom,'NELEM_SEND'    ,trim(vacal),nelem_send)
    call memory_deallo(memor_dom,'NPOIN_SEND'    ,trim(vacal),npoin_send)
    call memory_deallo(memor_dom,'NPOIN_RECV'    ,trim(vacal),npoin_recv)
    call memory_deallo(memor_dom,'NELEM_RECV'    ,trim(vacal),nelem_recv)

    do ipart = 1,npart 
       call geome_send(ipart) % deallo()
       call geome_recv(ipart) % deallo()
    end do
    
    mesh_new % comm % SIZE4          = SIZE4  
    mesh_new % comm % RANK4          = RANK4  
    mesh_new % comm % PAR_COMM_WORLD = int(PAR_COMM4,ip)

  end subroutine AMR_new_mesh
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-02
  !> @brief   Domain
  !> @details Recompute the domain
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_update_original_mesh(mesh_complete)
    
    type(mesh_type), intent(inout) :: mesh_complete

    call messages_live('UPDATE DOMAIN')
    !
    ! Copy new mesh
    !
    nelem = mesh_complete % nelem
    npoin = mesh_complete % npoin
    nboun = mesh_complete % nboun
    mnode = mesh_complete % mnode
    call PAR_MAX(mnode)
    !
    ! Deallocate old mesh
    !
    call memory_deallo(memor_dom,'LNODS'    ,'mesh_type_save_original_mesh',lnods    )
    call memory_deallo(memor_dom,'LTYPE'    ,'mesh_type_save_original_mesh',ltype    )
    call memory_deallo(memor_dom,'LEINV_LOC','mesh_type_save_original_mesh',leinv_loc)

    call memory_deallo(memor_dom,'COORD'    ,'mesh_type_save_original_mesh',coord    )
    call memory_deallo(memor_dom,'LNINV_LOC','mesh_type_save_original_mesh',lninv_loc)
    
    call memory_deallo(memor_dom,'LNODB'    ,'mesh_type_save_original_mesh',lnodb    )
    call memory_deallo(memor_dom,'LTYPB'    ,'mesh_type_save_original_mesh',ltypb    )
    call memory_deallo(memor_dom,'LELBO'    ,'mesh_type_save_original_mesh',lelbo    )
    call memory_deallo(memor_dom,'LBOCH'    ,'mesh_type_save_original_mesh',lboch    )
    call memory_deallo(memor_dom,'LBINV_LOC','mesh_type_save_original_mesh',lbinv_loc)
    !
    ! Copy new mesh from MESH_COMPLETE
    ! 
    call memory_copy(memor_dom,'MESH_COMPLETE % LNODS'    ,'mesh_type_save_original_mesh',mesh_complete % lnods    ,lnods    ,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_COMPLETE % LTYPE'    ,'mesh_type_save_original_mesh',mesh_complete % ltype    ,ltype    ,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_COMPLETE % LEINV_LOC','mesh_type_save_original_mesh',mesh_complete % leinv_loc,leinv_loc,'DO_NOT_DEALLOCATE')

    call memory_copy(memor_dom,'MESH_COMPLETE % COORD'    ,'mesh_type_save_original_mesh',mesh_complete % coord    ,coord    ,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_COMPLETE % LNINV_LOC','mesh_type_save_original_mesh',mesh_complete % lninv_loc,lninv_loc,'DO_NOT_DEALLOCATE')

    call memory_copy(memor_dom,'MESH_COMPLETE % LNODB'    ,'mesh_type_save_original_mesh',mesh_complete % lnodb    ,lnodb    ,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_COMPLETE % LTYPB'    ,'mesh_type_save_original_mesh',mesh_complete % ltypb    ,ltypb    ,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_COMPLETE % LELBO'    ,'mesh_type_save_original_mesh',mesh_complete % lelbo    ,lelbo    ,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_COMPLETE % LNOCH'    ,'mesh_type_save_original_mesh',mesh_complete % lboch    ,lboch    ,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_COMPLETE % LBINV_LOC','mesh_type_save_original_mesh',mesh_complete % lbinv_loc,lbinv_loc,'DO_NOT_DEALLOCATE')
    !
    ! Deallocate others
    !
    call domain_memory_deallocate('LNNOD')
    call domain_memory_deallocate('LNNOB')
    call domain_memory_deallocate('LGAUS')
    call domain_memory_deallocate('LPOTY')
    call domain_memory_deallocate('LNLEV')
    
  end subroutine AMR_update_original_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-02
  !> @brief   Domain
  !> @details Create boundary mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_boundary_mesh(mesh_complete)

    type(mesh_type),       intent(inout) :: mesh_complete
    integer(ip)                          :: iboun,ii,ipoin

    call messages_live('CREATE BOUNDARY MESH')

    if( mesh_complete % nelem > 0 ) then
       !
       ! Creat local boundary meshes
       !
       call meshes_list_boundary_elements(&
            mesh_complete % nelem,mesh_complete % npoin,mesh_complete % lnods,mesh_complete % ltype,&
            mesh_complete % nboun,mesh_complete % lelbo,mesh_complete % lnodb,mesh_complete % ltypb,&
            memor_dom,&
            LELBO_NAME='LELBO',&
            LNODB_NAME='LNODB',&
            LTYPB_NAME='LTYPB')
       !
       ! Remove internal faces
       !
       call AMR_remove_internal_boundaries(mesh_complete)
       !
       ! Compute LBOCH and LBINV_LOC
       !
       call memory_alloca(memor_dom,'LBOCH'    ,'mesh_type_save_original_mesh',mesh_complete % lboch    ,mesh_complete % nboun)
       call memory_alloca(memor_dom,'LBINV_LOC','mesh_type_save_original_mesh',mesh_complete % lbinv_loc,mesh_complete % nboun)
       do iboun = 1,mesh_complete % nboun
          mesh_complete % lboch(iboun) = BOFEM
       end do
       
    end if
    
    call par_global_numbering_elements_boundaries(mesh_complete % nboun,mesh_complete % lbinv_loc)

  end subroutine AMR_boundary_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-02
  !> @brief   Communicator
  !> @details Recompute communication arrays
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_update_communication_arrays(mesh_complete)
    
    type(mesh_type), intent(inout) :: mesh_complete
    
    call messages_live('UPDATE COMMUNICATION ARRAYS')
    call PAR_DEALLOCATE_COMMUNICATION_ARRAY(PAR_COMM_MY_CODE_ARRAY(1),par_memor)
    call PAR_COPY_COMMUNICATION_ARRAY(mesh_complete % comm,PAR_COMM_MY_CODE_ARRAY(1),par_memor)
    
    PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD =  PAR_COMM_MY_CODE
    PAR_COMM_MY_CODE_ARRAY(1) % RANK4          =  int(kfl_paral,4)
    commd                                      => PAR_COMM_MY_CODE_ARRAY(1)
    !
    ! Determine own nodes
    !
    call memory_copy(memor_dom,'BOUND_PERM','par_renumber_nodes',&
         PAR_COMM_MY_CODE_ARRAY(1) % bound_perm,&
         PAR_COMM_MY_CODE_ARRAY(1) % bound_invp,'DO_NOT_DEALLOCATE')
    
  end subroutine AMR_update_communication_arrays

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-02
  !> @brief   Domain
  !> @details Recompute the domain
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_domain(mesh_new)

    type(mesh_type),         intent(inout) :: mesh_new  
    integer(ip)                            :: ii,kk,ineig,lsize,dom_i
    integer(ip)                            :: kpoin,ipoin,ielem,pnode
    integer(ip)                            :: inode,iboun,ifiel
    integer(ip)                            :: color_target
    integer(ip)                            :: color_source
    integer(ip)                            :: PAR_CURRENT_RANK
    integer(ip)                            :: PAR_CURRENT_SIZE
    integer(ip), pointer                   :: ii1_source(:)
    integer(ip), pointer                   :: ii1_target(:)
    real(rp),    pointer                   :: xx1_source(:)
    real(rp),    pointer                   :: xx1_target(:)

    call messages_live('ARRAY INTERPOLATIONS')
    !
    ! Nullify pointers
    !
    nullify(ii1_source)
    nullify(ii1_target)
    nullify(xx1_source)
    nullify(xx1_target)

    !--------------------------------------------------------------------
    !
    ! Initialize couplings
    !
    !--------------------------------------------------------------------

    call cou_initialization(coupling_AMR_npoin)
    call cou_initialization(coupling_AMR_npoin_nearest)
    call cou_initialization(coupling_AMR_nelem)
    call cou_initialization(coupling_AMR_nboun)

    !--------------------------------------------------------------------
    !
    ! Create nodal coupling for element interplation
    !
    !--------------------------------------------------------------------

    coupling_AMR_npoin % itype                  =  ELEMENT_INTERPOLATION
    coupling_AMR_npoin % kfl_multi_source       =  0                          ! Enable multiple source nodes
    coupling_AMR_npoin % kfl_lost_wet_points    =  1                          ! Do not crash if wet points are lost
    coupling_AMR_npoin % number                 =  2010
    coupling_AMR_npoin % color_source           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    coupling_AMR_npoin % color_target           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    coupling_AMR_npoin % commd % PAR_COMM_WORLD =  commd % PAR_COMM_WORLD
    coupling_AMR_npoin % wet % npoin_wet        =  mesh_new % npoin
    coupling_AMR_npoin % target_entity          =  FLOATING_TARGET_ENTITY     ! Where coupling is eventuall applied
    color_target                                =  coupling_AMR_npoin % color_target
    color_source                                =  coupling_AMR_npoin % color_source
    !
    ! Define wet geometry: select all nodes
    !
    call memory_alloca(par_memor,'COUPLING % WET   % COORD_WET'  ,'par_interface_exchange',coupling_AMR_npoin % wet % coord_wet,ndime,coupling_AMR_npoin % wet % npoin_wet)
    call memory_alloca(par_memor,'COUPLING % COMMD % LRECV_PERM' ,'par_interface_exchange',coupling_AMR_npoin % commd % lrecv_perm,mesh_new % npoin)
    do ipoin = 1,mesh_new % npoin
       coupling_AMR_npoin % wet % coord_wet(:,ipoin)  = mesh_new % coord(:,ipoin)
       coupling_AMR_npoin % commd % lrecv_perm(ipoin) = ipoin
    end do
    !
    ! Coupling
    !
    call COU_INIT_INTERPOLATE_POINTS_VALUES(coupling_AMR_npoin,PAR_COMM_MY_CODE)
    call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling_AMR_npoin)

    !--------------------------------------------------------------------
    !
    ! Create nodal coupling for nearest node
    !
    !--------------------------------------------------------------------

    coupling_AMR_npoin_nearest % itype                  =  NEAREST_ELEMENT_NODE
    coupling_AMR_npoin_nearest % kfl_multi_source       =  0                          ! Enable multiple source nodes
    coupling_AMR_npoin_nearest % kfl_lost_wet_points    =  1                          ! Do not crash if wet points are lost
    coupling_AMR_npoin_nearest % number                 =  2013
    coupling_AMR_npoin_nearest % color_source           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    coupling_AMR_npoin_nearest % color_target           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    coupling_AMR_npoin_nearest % commd % PAR_COMM_WORLD =  commd % PAR_COMM_WORLD
    coupling_AMR_npoin_nearest % wet % npoin_wet        =  mesh_new % npoin
    coupling_AMR_npoin_nearest % target_entity          =  FLOATING_TARGET_ENTITY     ! Where coupling is eventuall applied
    color_target                                        =  coupling_AMR_npoin_nearest % color_target
    color_source                                        =  coupling_AMR_npoin_nearest % color_source
    !
    ! Define wet geometry: select all nodes
    !
    call memory_alloca(par_memor,'COUPLING % WET   % COORD_WET'  ,'par_interface_exchange',coupling_AMR_npoin_nearest % wet % coord_wet,ndime,coupling_AMR_npoin_nearest % wet % npoin_wet)
    call memory_alloca(par_memor,'COUPLING % COMMD % LRECV_PERM' ,'par_interface_exchange',coupling_AMR_npoin_nearest % commd % lrecv_perm,mesh_new % npoin)
    do ipoin = 1,mesh_new % npoin
       coupling_AMR_npoin_nearest % wet % coord_wet(:,ipoin)  = mesh_new % coord(:,ipoin)
       coupling_AMR_npoin_nearest % commd % lrecv_perm(ipoin) = ipoin
    end do
    !
    ! Coupling
    !
    call COU_INIT_INTERPOLATE_POINTS_VALUES(coupling_AMR_npoin_nearest,PAR_COMM_MY_CODE)
    call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling_AMR_npoin_nearest)

    !--------------------------------------------------------------------
    !
    ! Create element coupling
    !
    !--------------------------------------------------------------------

    coupling_AMR_nelem % itype                  =  ELEMENT_INTERPOLATION
    coupling_AMR_nelem % kfl_multi_source       =  0                          ! Enable multiple source nodes
    coupling_AMR_nelem % kfl_lost_wet_points    =  1                          ! Do not crash if wet points are lost
    coupling_AMR_nelem % number                 =  2011
    coupling_AMR_nelem % color_source           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    coupling_AMR_nelem % color_target           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    coupling_AMR_nelem % commd % PAR_COMM_WORLD =  commd % PAR_COMM_WORLD
    coupling_AMR_nelem % wet % npoin_wet        =  mesh_new % nelem
    coupling_AMR_nelem % target_entity          =  FLOATING_TARGET_ENTITY     ! Where coupling is eventuall applied
    color_target                                =  coupling_AMR_nelem % color_target
    color_source                                =  coupling_AMR_nelem % color_source
    coupling_AMR_nelem % kfl_source_value       =  VALUES_ON_ELEMENTS
    !
    ! Define wet geometry: select boundary nodes only
    !
    call memory_alloca(par_memor,'COUPLING % WET   % COORD_WET' ,'par_interface_exchange',coupling_AMR_nelem % wet % coord_wet,ndime,coupling_AMR_nelem % wet % npoin_wet)
    call memory_alloca(par_memor,'COUPLING % COMMD % LRECV_PERM','par_interface_exchange',coupling_AMR_nelem % commd % lrecv_perm,mesh_new % nelem)
    do ielem = 1,mesh_new % nelem
       coupling_AMR_nelem % wet % coord_wet(:,ielem) = 0.0_rp
       pnode = nnode(abs(mesh_new % ltype(ielem)))
       do inode = 1,pnode
          ipoin = mesh_new % lnods(inode,ielem)
          coupling_AMR_nelem % wet % coord_wet(:,ielem) = coupling_AMR_nelem % wet % coord_wet(:,ielem) + mesh_new % coord(:,ipoin)
       end do
       coupling_AMR_nelem % wet % coord_wet(:,ielem) = coupling_AMR_nelem % wet % coord_wet(:,ielem) / real(pnode,rp)
       coupling_AMR_nelem % commd % lrecv_perm(ielem) = ielem
    end do
    !
    ! Coupling
    !    
    call COU_INIT_INTERPOLATE_POINTS_VALUES(coupling_AMR_nelem) 
    call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling_AMR_nelem)

    !--------------------------------------------------------------------
    !
    ! Create boundary coupling
    !
    !--------------------------------------------------------------------

    coupling_AMR_nboun % itype                  =  BOUNDARY_INTERPOLATION
    coupling_AMR_nboun % kfl_multi_source       =  0                          ! Enable multiple source nodes
    coupling_AMR_nboun % kfl_lost_wet_points    =  1                          ! Do not crash if wet points are lost
    coupling_AMR_nboun % number                 =  2012
    coupling_AMR_nboun % color_source           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    coupling_AMR_nboun % color_target           =  par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    coupling_AMR_nboun % commd % PAR_COMM_WORLD =  commd % PAR_COMM_WORLD
    coupling_AMR_nboun % wet % npoin_wet        =  mesh_new % nboun
    coupling_AMR_nboun % target_entity          =  FLOATING_TARGET_ENTITY     ! Where coupling is eventuall applied
    color_target                                =  coupling_AMR_nboun % color_target
    color_source                                =  coupling_AMR_nboun % color_source
    coupling_AMR_nboun % kfl_source_value       =  VALUES_ON_BOUNDARIES
    !
    ! Define wet geometry: select boundary nodes only
    !
    call memory_alloca(par_memor,'COUPLING % WET % COORD_WET' ,'par_interface_exchange',coupling_AMR_nboun % wet % coord_wet,ndime,coupling_AMR_nboun % wet % npoin_wet)
    call memory_alloca(par_memor,'COUPLING % WET % LPOIN_WET' ,'par_interface_exchange',coupling_AMR_nboun % wet % lpoin_wet,coupling_AMR_nboun % wet % npoin_wet)
    do iboun = 1,mesh_new % nboun
       coupling_AMR_nboun % wet % coord_wet(:,iboun) = 0.0_rp
       pnode = nnode(abs(mesh_new % ltypb(iboun)))
       do inode = 1,pnode
          ipoin = mesh_new % lnodb(inode,iboun)
          coupling_AMR_nboun % wet % coord_wet(:,iboun) = coupling_AMR_nboun % wet % coord_wet(:,iboun) + mesh_new % coord(:,ipoin)
       end do
       coupling_AMR_nboun % wet % coord_wet(:,iboun) = coupling_AMR_nboun % wet % coord_wet(:,iboun) / real(pnode,rp)
       coupling_AMR_nboun % wet % lpoin_wet(iboun)   = iboun
    end do
    !
    ! Construct KD-Tree
    !
    call kdtree_construct(&
         nboun,npoin,lnodb,ltypb,coord,coupling_AMR_nboun % geome % kdtree)
    !
    ! Coupling
    !    
    call COU_INIT_INTERPOLATE_POINTS_VALUES(coupling_AMR_nboun)
    call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling_AMR_nboun)

    npoin_new = mesh_new % npoin
    nboun_new = mesh_new % nboun
    nelem_new = mesh_new % nelem

    !--------------------------------------------------------------------
    !
    ! Interpolate domain arrays 
    !
    !--------------------------------------------------------------------
    !
    ! Node boujndary conditions
    !
    !call domain_memory_deallocate('KFL_CODNO')
    !kfl_icodn = 0

    call AMR_interpolate(lnoch    ,'NPOIN',1_ip,memor_dom,'LNOCH')
    call AMR_interpolate(lmast    ,'NPOIN',1_ip,memor_dom,'LMAST')
    call AMR_interpolate(lnset    ,'NPOIN',1_ip,memor_dom,'LNSET')
    call AMR_interpolate(lgrou_dom,'NPOIN',1_ip,memor_dom,'LGROU_DOM')
    call AMR_interpolate(kfl_codno,'NPOIN',2_ip,memor_dom,'KFL_CODNO')

    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
          if( kfl_field(4,ifiel) == 1 ) then
             call AMR_interpolate(xfiel(ifiel) % a,'NPOIN',2_ip,memor_dom,VARIABLE_NAME='XFIEL % A')
          end if
       end if
    end do

    call AMR_interpolate(lelch,    'NELEM',1_ip,memor_dom,'LELCH')
    call AMR_interpolate(lesub,    'NELEM',1_ip,memor_dom,'LESUB')
    call AMR_interpolate(lmate,    'NELEM',1_ip,memor_dom,'LMATE')
    call AMR_interpolate(leset,    'NELEM',1_ip,memor_dom,'LESET')
    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NELEM_TYPE ) then
          if( kfl_field(4,ifiel) == 1 ) then
             call AMR_interpolate(xfiel(ifiel) % a,'NELEM',2_ip,memor_dom,VARIABLE_NAME='XFIEL % A')
          end if
       end if
    end do

    call AMR_interpolate(kfl_codbo,'NBOUN',2_ip,memor_dom,'KFL_CODNO')
    call AMR_interpolate(lbset,    'NBOUN',1_ip,memor_dom,'LBSET')
    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
          if( kfl_field(4,ifiel) == 1 ) then
             call AMR_interpolate(xfiel(ifiel) % a,'NBOUN',2_ip,memor_dom,VARIABLE_NAME='XFIEL % A')
          end if
       end if
    end do
    
    nelem = mesh_new % nelem
    npoin = mesh_new % npoin
    nboun = mesh_new % nboun

    call par_global_variables_arrays()
    call par_mesh_dimensions()

  end subroutine AMR_domain

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-09
  !> @brief   Remove internal boundaries
  !> @details Identify common faces with neighbors and remove them
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_remove_internal_boundaries(mesh)

    type(mesh_type),     intent(inout) :: mesh
    integer(ip),         pointer       :: nboun_send(:)
    type(i2p),           pointer       :: lboun_send(:)
    integer(ip),         pointer       :: nboun_recv(:)
    type(i2p),           pointer       :: lboun_recv(:)
    integer(ip),         pointer       :: nsubd_npoin(:)
    type(i1p),           pointer       :: lsubd_npoin(:)
    integer(ip),         pointer       :: lboun_tot(:,:)
    integer(ip),         pointer       :: knodb(:)
    integer(ip)                        :: ii,ipoin,iboun,inodb,pnodb
    integer(ip)                        :: ineig,dom_i,kneig,kboun
    integer(ip)                        :: lnodb_loc(mnodb)
    integer(4)                         :: PAR_COMM_AMR4
    logical(lg)                        :: same_boundary
    integer(ip),         pointer       :: lnodb_cpy(:,:)   
    integer(ip),         pointer       :: ltypb_cpy(:)    
    integer(ip),         pointer       :: lelbo_cpy(:)  

    if( mesh % comm % nneig > 0 ) then

       nullify(nboun_send)
       nullify(lboun_send)
       nullify(nboun_recv)
       nullify(lboun_recv)
       nullify(nsubd_npoin)
       nullify(lsubd_npoin)
       nullify(lboun_tot)
       nullify(knodb)

       nullify(lnodb_cpy)
       nullify(ltypb_cpy)
       nullify(lelbo_cpy)

       PAR_COMM_AMR4 = mesh % comm % PAR_COMM_WORLD

       call memory_alloca(memor_dom,'NBOUN_SEND' ,trim(vacal),nboun_send ,mesh % comm % nneig)
       call memory_alloca(memor_dom,'LBOUN_SEND' ,trim(vacal),lboun_send ,mesh % comm % nneig)
       call memory_alloca(memor_dom,'LBOUN_TOT'  ,trim(vacal),lboun_tot  ,mesh % comm % nneig,mesh % nboun)
       call memory_alloca(memor_dom,'NSUBD_NPOIN',trim(vacal),nsubd_npoin,mesh % npoin)
       call memory_alloca(memor_dom,'LSUBD_NPOIN',trim(vacal),lsubd_npoin,mesh % npoin)
       do ipoin = 1,mesh % npoin
          call memory_alloca(memor_dom,'LSUBD_NPOIN % L',trim(vacal),lsubd_npoin(ipoin)%l,mesh % comm % nneig)
       end do
       call memory_alloca(memor_dom,'KNODB' ,trim(vacal),knodb,mesh % comm % nneig)
       !
       ! List of node neighbors LSUBD_NPOIN(1:NPOIN) % L(:) 
       !
       do ineig = 1,mesh % comm % nneig
          do ii = mesh % comm % bound_size(ineig),mesh % comm % bound_size(ineig+1)-1
             ipoin                                      = mesh % comm % bound_perm(ii)
             nsubd_npoin(ipoin)                         = nsubd_npoin(ipoin) + 1
             lsubd_npoin(ipoin) % l(nsubd_npoin(ipoin)) = ineig
          end do
       end do
       !
       ! NBOUN_SEND(INEIG): Number of boundaries to send
       !
       do iboun = 1,mesh % nboun
          knodb = 0
          pnodb = nnode(abs(mesh % ltypb(iboun)))
          do inodb = 1,pnodb
             ipoin = mesh % lnodb(inodb,iboun)
             do ii = 1,nsubd_npoin(ipoin)
                ineig        = lsubd_npoin(ipoin) % l(ii)
                knodb(ineig) = knodb(ineig) + 1
             end do
          end do
          kneig = 0
          do ineig = 1,mesh % comm % nneig
             if( knodb(ineig) == pnodb ) then
                kneig                  = kneig + 1
                nboun_send(ineig)      = nboun_send(ineig) + 1
                lboun_tot(kneig,iboun) = ineig
             end if
          end do
       end do
       call memory_deallo(memor_dom,'NSUBD_NPOIN',trim(vacal),nsubd_npoin)
       call memory_deallo(memor_dom,'LSUBD_NPOIN',trim(vacal),lsubd_npoin)
       !
       ! LBOUN_SEND(INEIG) % L(:,:): List of boundaries to send
       !
       do ineig = 1,mesh % comm % nneig
          call memory_alloca(memor_dom,'LBOUN_SEND % L',trim(vacal),lboun_send(ineig) % l,mesh % mnodb,nboun_send(ineig))
          nboun_send(ineig) = 0
       end do

       do iboun = 1,mesh % nboun
          loop_kneig: do kneig = 1,mesh % comm % nneig
             ineig = lboun_tot(kneig,iboun)
             if( ineig == 0 ) then
                exit loop_kneig
             else
                pnodb             = nnode(abs(mesh % ltypb(iboun)))
                nboun_send(ineig) = nboun_send(ineig) + 1
                do inodb = 1,pnodb
                   ipoin = mesh % lnodb(inodb,iboun)
                   lboun_send(ineig) % l(inodb,nboun_send(ineig)) = mesh % lninv_loc(ipoin)
                end do
                call maths_heap_sort(2_ip,pnodb,lboun_send(ineig) % l(:,nboun_send(ineig)))
             end if
          end do loop_kneig
       end do
       !
       ! Send/Receive list
       !
       call memory_alloca(memor_dom,'NBOUN_RECV' ,trim(vacal),nboun_recv ,mesh % comm % nneig)
       call memory_alloca(memor_dom,'LBOUN_RECV' ,trim(vacal),lboun_recv ,mesh % comm % nneig)

       do ineig = 1,mesh % comm % nneig
          dom_i = mesh % comm % neights(ineig)
          call PAR_SEND_RECEIVE(nboun_send(ineig),nboun_recv(ineig),'IN MY CODE',dom_i,PAR_COMM_IN4=PAR_COMM_AMR4)
          call memory_alloca(memor_dom,'LBOUN_RECV % L',trim(vacal),lboun_recv(ineig) % l,mesh % mnodb,nboun_recv(ineig))
          call PAR_SEND_RECEIVE(lboun_send(ineig) %l,lboun_recv(ineig) %l,'IN MY CODE',dom_i,PAR_COMM_IN4=PAR_COMM_AMR4)
       end do
       !
       ! Identify internal boundaries
       !
       do iboun = 1,mesh % nboun
          pnodb              = nnode(abs(mesh % ltypb(iboun)))
          lnodb_loc(1:pnodb) = mesh % lninv_loc(mesh % lnodb(1:pnodb,iboun))
          lboun_tot(1,iboun) = 1
          call maths_heap_sort(2_ip,pnodb,lnodb_loc)
          loop_ineig: do ineig = 1,mesh % comm % nneig
             do kboun = 1,nboun_recv(ineig)
                same_boundary = .true.
                loop_inodb: do inodb = 1,pnodb
                   if( lnodb_loc(inodb) /= lboun_recv(ineig) %l(inodb,kboun) ) then
                      same_boundary = .false.
                      exit loop_inodb
                   end if
                end do loop_inodb
                if( same_boundary ) then
                   lboun_tot(1,iboun) = 0
                   exit loop_ineig
                end if
             end do
          end do loop_ineig
       end do
       !
       ! Deallocate
       !
       call memory_deallo(memor_dom,'NBOUN_SEND',trim(vacal),nboun_send )
       call memory_deallo(memor_dom,'LBOUN_SEND',trim(vacal),lboun_send )
       call memory_deallo(memor_dom,'KNODB'     ,trim(vacal),knodb      )
       !
       ! Resize boundary arrays
       !
       kboun = 0
       do iboun = 1,mesh % nboun
          kboun = kboun + lboun_tot(1,iboun)       
       end do
       call memory_copy  (memor_dom,'LNODB_CPY'           ,trim(vacal),mesh % lnodb    ,lnodb_cpy)
       call memory_copy  (memor_dom,'LTYPB_CPY'           ,trim(vacal),mesh % ltypb    ,ltypb_cpy)
       call memory_copy  (memor_dom,'LELBO_CPY'           ,trim(vacal),mesh % lelbo    ,lelbo_cpy)
       !call memory_alloca(memor_dom,'MESH % LNODB'    ,trim(vacal),mesh % lnodb    ,mesh % mnodb,kboun,'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'MESH % LNODB'    ,trim(vacal),mesh % lnodb    ,mnodb,kboun,'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'MESH % LTYPB'    ,trim(vacal),mesh % ltypb    ,kboun                 ,'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,'MESH % LELBO'    ,trim(vacal),mesh % lelbo    ,kboun                 ,'DO_NOT_INITIALIZE')

       kboun = 0
       do iboun = 1,mesh % nboun
          if( lboun_tot(1,iboun) == 1 ) then
             kboun                           = kboun + 1
             pnodb                           = nnode(abs(ltypb_cpy(iboun)))
             mesh % lnodb(1:pnodb,kboun) = lnodb_cpy(1:pnodb,iboun)
             mesh % ltypb(kboun)         = ltypb_cpy(iboun)
             mesh % lelbo(kboun)         = lelbo_cpy(iboun)
          end if
       end do
       mesh % nboun = kboun

       call memory_deallo(memor_dom,'LBOUN_TOT' ,trim(vacal),lboun_tot)
       call memory_deallo(memor_dom,'LNODB_CPY' ,trim(vacal),lnodb_cpy)
       call memory_deallo(memor_dom,'LTYPB_CPY' ,trim(vacal),ltypb_cpy)
       call memory_deallo(memor_dom,'LELBO_CPY' ,trim(vacal),lelbo_cpy)

    end if

  end subroutine AMR_remove_internal_boundaries

end module mod_AMR
!> @}
