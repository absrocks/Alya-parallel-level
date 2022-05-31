!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> Manage the memory coupling
!> @{
!> @file    mod_coupling_memory.f90
!> @author  houzeaux
!> @date    2018-03-27
!> @brief   Coupling memory
!> @details Manages the memory of the coupling structure
!>          
!-----------------------------------------------------------------------

module mod_coupling_memory

  use def_kintyp, only : ip,rp,lg
  use def_kintyp, only : spmat
  use def_coupli, only : typ_color_coupling
  use def_coupli, only : typ_coupling_geometry
  use def_coupli, only : typ_coupling_wet
  use mod_parall, only : PAR_INITIALIZE_COMMUNICATION_ARRAY
  use mod_parall, only : PAR_DEALLOCATE_COMMUNICATION_ARRAY
  use def_master, only : AT_BEGINNING
  use def_coupli, only : memor_cou
  use def_coupli, only : RELAXATION_SCHEME
  use def_coupli, only : VALUES_ON_NODES
  use def_coupli, only : UNKNOWN
  use mod_kdtree, only : kdtree_deallocate
  use mod_kdtree, only : kdtree_initialize
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  implicit none
  private

  interface cou_initialization
     module procedure cou_initialization_array,&
          &           cou_initialization_scalar 
  end interface cou_initialization
 
  interface cou_deallocate
     module procedure cou_deallocate_array,&
          &           cou_deallocate_scalar 
  end interface cou_deallocate
 
  public :: cou_initialization
  public :: cou_deallocate
  public :: COU_DEALLOCATE_GEOMETRY   
  public :: COU_DEALLOCATE_TRANSMISSION          
  public :: COU_DEALLOCATE_SINGLE_COUPLING
  public :: COU_INITIALIZATION_SINGLE_COUPLING
  
contains
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    12/03/2018
  !> @brief   Initialize coupling  
  !> @details Initialize coupling variables of all couplings
  !>
  !----------------------------------------------------------------------
 
  subroutine cou_initialization_array(coupling)
 
    type(typ_color_coupling), intent(inout), pointer :: coupling(:)
    integer(ip) :: icoup
    
    do icoup = 1,size(coupling)
       call COU_INITIALIZATION_SINGLE_COUPLING(coupling(icoup))
    end do

  end subroutine cou_initialization_array

  subroutine cou_initialization_scalar(coupling)
 
    type(typ_color_coupling), intent(inout) :: coupling
    
    call COU_INITIALIZATION_SINGLE_COUPLING(coupling)

  end subroutine cou_initialization_scalar

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    12/03/2018
  !> @brief   Initialize coupling  
  !> @details Deallocate coupling structure
  !>
  !----------------------------------------------------------------------
 
  subroutine cou_deallocate_array(coupling)
 
    type(typ_color_coupling), intent(inout), pointer :: coupling(:)
    integer(ip) :: icoup
    
    do icoup = 1,size(coupling)
       call COU_DEALLOCATE_SINGLE_COUPLING(coupling(icoup))
    end do

  end subroutine cou_deallocate_array

  subroutine cou_deallocate_scalar(coupling)
 
    type(typ_color_coupling), intent(inout) :: coupling
    
    call COU_DEALLOCATE_SINGLE_COUPLING(coupling)

  end subroutine cou_deallocate_scalar

  !---------------------------------------------------------------------- 
  !
  !> @author  Matias Rivero
  !> @date    16/12/2014
  !> @brief   Deallocate coupling structure (geome)
  !> @details Deallocate coupling structure (geome)
  !>
  !----------------------------------------------------------------------

  subroutine COU_DEALLOCATE_WET(wet)

    type(typ_coupling_wet), intent(inout) :: wet
    integer(ip)                           :: iboun_wet

    wet % number_wet_points = 0
    wet % npoin_wet         = 0
    wet % nboun_wet         = 0
    
    call memory_deallo(memor_cou,'COUPLING % WET % VMASB_WET'     ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % vmasb_wet  )
    call memory_deallo(memor_cou,'COUPLING % WET % LBOUN_WET'     ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % lboun_wet  )
    call memory_deallo(memor_cou,'COUPLING % WET % LPOIN_WET'     ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % lpoin_wet  )
    call memory_deallo(memor_cou,'COUPLING % WET % COORD_WET'     ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % coord_wet  )
    call memory_deallo(memor_cou,'COUPLING % WET % WEIGHT_WET'    ,'DEALLOCATE_COUPLING_STRUCTURE_WET',wet % weight_wet )
    call memory_deallo(memor_cou,'COUPLING % WET % WEIGHT_WET_IMP','DEALLOCATE_COUPLING_STRUCTURE_WET',wet % weight_wet_imp )
    call memory_deallo(memor_cou,'COUPLING % WET % MASS_RATIO_WET','DEALLOCATE_COUPLING_STRUCTURE_WET',wet % mass_ratio_wet )

    if( associated(wet % proje_target) ) then
       do iboun_wet = 1,size(wet % proje_target)
          if( associated(wet % proje_target(iboun_wet) % permu) ) deallocate(wet % proje_target(iboun_wet) % permu)
          if( associated(wet % proje_target(iboun_wet) % shapb) ) deallocate(wet % proje_target(iboun_wet) % shapb)
          if( associated(wet % proje_target(iboun_wet) % gbsur) ) deallocate(wet % proje_target(iboun_wet) % gbsur)
       end do
       deallocate(wet % proje_target)
    end if

  end subroutine COU_DEALLOCATE_WET

  !---------------------------------------------------------------------- 
  !
  !> @author  Matias Rivero
  !> @date    16/12/2014
  !> @brief   Deallocate coupling structure (geome)
  !> @details Deallocate coupling structure (geome)
  !>
  !----------------------------------------------------------------------

  subroutine COU_DEALLOCATE_GEOMETRY(geome)

    type(typ_coupling_geometry), intent(inout) :: geome

    geome % ndime     = 0

    geome % nelem_source = 0
    geome % nboun_source = 0
    geome % npoin_source = 0

    call memory_deallo(memor_cou,'COUPLING % GEOME % LELEM_SOURCE', 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % lelem_source )
    call memory_deallo(memor_cou,'COUPLING % GEOME % LBOUN_SOURCE', 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % lboun_source )
    call memory_deallo(memor_cou,'COUPLING % GEOME % LPOIN_SOURCE', 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % lpoin_source )
    call memory_deallo(memor_cou,'COUPLING % GEOME % STATUS'      , 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % status       )
    call memory_deallo(memor_cou,'COUPLING % GEOME % SHAPF'       , 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % shapf        )    
    call memory_deallo(memor_cou,'COUPLING % GEOME % VMASB'       , 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % vmasb        )
    call memory_deallo(memor_cou,'COUPLING % GEOME % SCHED_PERM'  , 'DEALLOCATE_COUPLING_STRUCTURE_GEOME',geome % sched_perm   )

    call kdtree_deallocate(geome % kdtree)

  end subroutine COU_DEALLOCATE_GEOMETRY

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize one single coupling 
  !> @details Initialize one single coupling 
  !>
  !----------------------------------------------------------------------
 
  subroutine COU_INITIALIZATION_SINGLE_COUPLING(coupling,WET,INPUT)

    use def_coupli,      only : NODE_TARGET_ENTITY

    type(typ_color_coupling), intent(inout)        :: coupling
    logical(lg),              intent(in), optional :: WET
    logical(lg),              intent(in), optional :: INPUT
    logical(lg)                                    :: ifwet
    logical(lg)                                    :: ifinput

    ifwet = .true.
    if( present(WET) ) ifwet = WET

    ifinput = .true.
    if( present(INPUT) ) ifinput = INPUT

    if( ifinput ) then
       coupling % number                   =  0 
       coupling % itype                    =  0 
       coupling % code_source              =  1 
       coupling % code_target              =  1 
       coupling % zone_source              =  0 
       coupling % zone_target              =  0 
       coupling % color_source             =  0 
       coupling % color_target             =  0 
       coupling % module_source            =  0 
       coupling % module_target            =  0 
       coupling % subdomain_source         =  0 
       coupling % subdomain_target         =  0 
       coupling % where_type               =  0 
       coupling % where_number             =  0
       coupling % where_type_source        =  0 
       coupling % where_number_source      =  0
       coupling % target_entity            =  NODE_TARGET_ENTITY 
       coupling % kfl_toda_costa           =  1 
       coupling % kfl_par_transm           =  0 
       coupling % kfl_check_exha           =  0 
       coupling % what                     =  UNKNOWN 
       coupling % scheme                   =  RELAXATION_SCHEME     
       coupling % itera                    =  0 
       coupling % conservation             =  0 
       coupling % overlap                  =  0            
       coupling % ngaus                    =  0
       coupling % kfl_symmetry             =  0
       coupling % kfl_source_value         =  VALUES_ON_NODES
       coupling % kfl_multi_source         =  0
       coupling % kfl_lost_wet_points      =  0
       coupling % kind                     =  0      
       coupling % mirror_coupling          =  0      
       coupling % task_compute_and_send    = -100_ip 
       coupling % when_compute_and_send    = -100_ip 
       coupling % task_recv_and_assemble   = -100_ip 
       coupling % when_recv_and_assemble   = -100_ip 
       coupling % frequ_send               =  0      
       coupling % frequ_recv               =  0      
       coupling % when_update_wet_geometry =  AT_BEGINNING
       coupling % when_update_coupling     =  AT_BEGINNING
       coupling % temporal_predictor       =  0      
       coupling % temporal_predictor_order = -1_ip      
       coupling % variable                 =  '     '
       coupling % relax                    =  1.0_rp
       coupling % resid                    =  0.0_rp
       coupling % cputim                   =  0.0_rp
       coupling % aitken                   =  1.0_rp
       coupling % ranku_iqnls              = -1_ip
       coupling % history_iqnls            =  1_ip
       coupling % scaling_iqnls            = -1.0_rp
       coupling % efilter_iqnls            = -1.0_rp
    end if
              
    nullify(coupling % values)           
    nullify(coupling % values_frequ)     
    nullify(coupling % values_predicted)   
    nullify(coupling % values_converged) 
    nullify(coupling % jacobian_inverse) 
    nullify(coupling % dincr_predicted)      
    nullify(coupling % residues_iqnls)    
    nullify(coupling % relaxed_iqnls)     
    nullify(coupling % unrelaxed_iqnls)   
    nullify(coupling % valincr_iqnls)     
    nullify(coupling % residincr_iqnls)
    nullify(coupling % valincr_history_iqnls)
    nullify(coupling % residincr_history_iqnls)
    nullify(coupling % V_current_history_iqnls)
    nullify(coupling % W_current_history_iqnls)
    nullify(coupling % history_tracking_iqnls)

    call COU_INITIALIZATION_GEOMETRY       (coupling % geome)
    if( ifwet ) call COU_INITIALIZATION_WET(coupling % wet)
    call PAR_INITIALIZE_COMMUNICATION_ARRAY(coupling % commd,PAR_COMM_OPT=.false.)
    call COU_INITIALIZATION_TRANSMISSION   (coupling % ltransmat_target)

  end subroutine COU_INITIALIZATION_SINGLE_COUPLING

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize coupling 
  !> @details Initialize coupling transmission matrices 
  !>
  !----------------------------------------------------------------------

  subroutine COU_INITIALIZATION_TRANSMISSION(ltransmat)
    
    type(spmat), pointer, intent(inout) :: ltransmat(:) ! Transmission matrix
    integer(ip)                         :: ii

    nullify(ltransmat)
    !if( associated(ltransmat) ) then
    !   do ii = 1,size(ltransmat)
    !      ltransmat(ii) % ndof  = 0
    !      ltransmat(ii) % nrows = 0
    !      ltransmat(ii) % ncols = 0
    !      
    !      nullify(ltransmat(ii) % iA)
    !      nullify(ltransmat(ii) % jA)
    !      nullify(ltransmat(ii) % vA)
    !   end do
    !end if
    
  end subroutine COU_INITIALIZATION_TRANSMISSION
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize coupling 
  !> @details Initialize coupling geometry
  !>
  !----------------------------------------------------------------------

  subroutine COU_INITIALIZATION_WET(wet)

    use def_coupli,      only : NODE_WET_POINT
    type(typ_coupling_wet), intent(inout) :: wet
          
     wet % number_wet_points = 0  
     wet % npoin_wet         = 0   
     wet % nboun_wet         = 0
     wet % point_type        = NODE_WET_POINT

     
     nullify(wet % lboun_wet)               
     nullify(wet % lpoin_wet)               
     nullify(wet % coord_wet)             
     nullify(wet % weight_wet)              
     nullify(wet % weight_wet_imp)              
     nullify(wet % vmasb_wet)               
     nullify(wet % mass_ratio_wet)          
     nullify(wet % proje_target)

  end subroutine COU_INITIALIZATION_WET
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize coupling 
  !> @details Initialize coupling geometry
  !>
  !----------------------------------------------------------------------

  subroutine COU_INITIALIZATION_GEOMETRY(geome)

    use def_coupli, only : typ_coupling_geometry
    type(typ_coupling_geometry), intent(inout) :: geome

     geome % ndime             = 0     
     geome % nelem_source      = 0      
     geome % nboun_source      = 0         
     geome % npoin_source      = 0 
     
     nullify(geome % lelem_source)            
     nullify(geome % lboun_source)            
     nullify(geome % lpoin_source)            
     nullify(geome % status)                  
     nullify(geome % shapf)               
     nullify(geome % vmasb)
     nullify(geome % sched_perm) 

     call kdtree_initialize(geome % kdtree) 

  end subroutine COU_INITIALIZATION_GEOMETRY

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize coupling 
  !> @details Initialize coupling transmission matrices 
  !>
  !----------------------------------------------------------------------

  subroutine COU_DEALLOCATE_TRANSMISSION(ltransmat)

    type(spmat), pointer, intent(inout) :: ltransmat(:)     ! Transmission matrix

    if( associated(ltransmat) ) then
       ltransmat(:) % ndof  = 0
       ltransmat(:) % nrows = 0
       ltransmat(:) % ncols = 0  
       call memory_deallo(memor_cou,'LTRANSMAT_TARGET','COU_DEALLOCATE_SINGLE_COUPLING',ltransmat)
    end if
    
  end subroutine COU_DEALLOCATE_TRANSMISSION
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/04/2014
  !> @brief   Initialize one single coupling 
  !> @details Initialize one single coupling 
  !>
  !----------------------------------------------------------------------
 
  subroutine COU_DEALLOCATE_SINGLE_COUPLING(coupling,WET)

    type(typ_color_coupling), intent(inout)        :: coupling
    logical(lg),              intent(in), optional :: WET
    logical(lg)                                    :: ifwet

    ifwet = .true.
    if( present(WET) ) ifwet = WET
    
    call memory_deallo(memor_cou,'coupling % values          ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % values          )
    call memory_deallo(memor_cou,'coupling % values_frequ    ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % values_frequ    )     
    call memory_deallo(memor_cou,'coupling % values_predicted','COU_DEALLOCATE_SINGLE_COUPLING',coupling % values_predicted)   
    call memory_deallo(memor_cou,'coupling % values_converged','COU_DEALLOCATE_SINGLE_COUPLING',coupling % values_converged) 
    call memory_deallo(memor_cou,'coupling % jacobian_inverse','COU_DEALLOCATE_SINGLE_COUPLING',coupling % jacobian_inverse) 
    call memory_deallo(memor_cou,'coupling % dincr_predicted ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % dincr_predicted )      
    call memory_deallo(memor_cou,'coupling % residues_iqnls  ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % residues_iqnls  )    
    call memory_deallo(memor_cou,'coupling % relaxed_iqnls   ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % relaxed_iqnls   )     
    call memory_deallo(memor_cou,'coupling % unrelaxed_iqnls ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % unrelaxed_iqnls )   
    call memory_deallo(memor_cou,'coupling % valincr_iqnls   ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % valincr_iqnls   )     
    call memory_deallo(memor_cou,'coupling % residincr_iqnls ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % residincr_iqnls )    
    call memory_deallo(memor_cou,'coupling % valincr_history_iqnls ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % valincr_history_iqnls )    
    call memory_deallo(memor_cou,'coupling % residincr_history_iqnls ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % residincr_history_iqnls )    
    call memory_deallo(memor_cou,'coupling % V_current_history_iqnls ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % V_current_history_iqnls )    
    call memory_deallo(memor_cou,'coupling % W_current_history_iqnls ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % W_current_history_iqnls )    
    call memory_deallo(memor_cou,'coupling % history_tracking_iqnls ','COU_DEALLOCATE_SINGLE_COUPLING',coupling % history_tracking_iqnls )    
    
    call COU_DEALLOCATE_GEOMETRY           (coupling % geome)
    if( ifwet ) call COU_DEALLOCATE_WET    (coupling % wet)
    call PAR_DEALLOCATE_COMMUNICATION_ARRAY(coupling % commd,memor_cou,PAR_COMM_OPT=.false.,COMM_NAME='COUPLING % COMMD',INITIALIZE=.true.)
    call COU_DEALLOCATE_TRANSMISSION       (coupling % ltransmat_target)   
    
  end subroutine COU_DEALLOCATE_SINGLE_COUPLING

end module mod_coupling_memory
!> @}
