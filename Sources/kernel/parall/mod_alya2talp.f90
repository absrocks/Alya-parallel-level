!-----------------------------------------------------------------------
!> @addtogroup Parall
!> Brige to TALP
!> @{
!> @file    mod_alya2talp.f90
!> @author  houzeaux
!> @date    2019-01-07
!> @brief   Bridge to TALP
!> @details Interfaces with TALP
!>          Data for Dynamic Allocation campaign:
!>          DLB_MONITOR_DA
!>          DLB_HANDLE_DA
!>          Data for global timings:
!>          DLB_MONITOR_GLOBAL
!>          DLB_HANDLE_GLOBAL
!>          Data for modules:
!>          DLB_MONITOR_MODULE
!>          DLB_HANDLE_MODULE
!>
!-----------------------------------------------------------------------

module mod_alya2talp

  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : mmodu,modul
  use def_master,         only : namod
  use def_master,         only : kfl_modul
  use def_master,         only : intost
  use def_master,         only : zeror
  use def_master,         only : npart
  use mod_messages,       only : messages_live
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_AVERAGE

  use, intrinsic :: ISO_C_BINDING

  implicit none
     
#if defined ALYA_TALP
  include 'dlbf_talp.h'

  type talp_module
     type(dlb_monitor_t), pointer :: dlb_monitor
     type(c_ptr)                  :: dlb_handle
  end type talp_module
     
  type(dlb_monitor_t), pointer :: dlb_monitor_da
  type(c_ptr)                  :: dlb_handle_da

  type(dlb_monitor_t), pointer :: dlb_monitor_global
  type(c_ptr)                  :: dlb_handle_global

  type(dlb_monitor_t), pointer :: dlb_monitor_timing
  type(c_ptr)                  :: dlb_handle_timing

  type(talp_module)            :: dlb_module(0:mmodu)

#endif  

  private

  public :: alya2talp_register
  public :: alya2talp_MonitoringRegionStart
  public :: alya2talp_MonitoringRegionStop
  public :: alya2talp_MonitoringRegionreset
  public :: alya2talp_elapsed_time
  public :: alya2talp_parallel_efficiency
  public :: alya2talp_initialization

#if defined ALYA_TALP
  public :: dlb_module
  public :: dlb_monitor_da
  public :: dlb_handle_da
#endif  

contains

  subroutine alya2talp_initialization()

    integer(ip) :: imodu

#if defined ALYA_TALP
    nullify(dlb_monitor_da)
    nullify(dlb_monitor_timing)
    nullify(dlb_monitor_global)
    do imodu = 0,mmodu
       nullify(dlb_module(imodu) % dlb_monitor)
    end do
#endif  

  end subroutine alya2talp_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2019-07-19
  !> @brief   Initialization
  !> @details Initialization of TALP whenevr TALP is used for load balance
  !>          or TALP_Barrier is used
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2talp_register()

    integer(ip) :: imodu

#if defined ALYA_TALP
    dlb_handle_da     = DLB_MonitoringRegionRegister(c_char_"region dynamic allocation"//C_NULL_CHAR)
    dlb_handle_timing = DLB_MonitoringRegionRegister(c_char_"region timing"//C_NULL_CHAR)
    dlb_handle_global = DLB_MonitoringRegionRegister(c_char_"region global"//C_NULL_CHAR)
    do imodu = 0,mmodu
       dlb_module(imodu) % dlb_handle = DLB_MonitoringRegionRegister(c_char_"region module "//trim(intost(imodu))//C_NULL_CHAR)
    end do

    if (.not. c_associated(dlb_handle_da))     call runend('MOD_ALYA2TALP: ERROR WHEN ASSOCIATING HANDLE_DA')
    if (.not. c_associated(dlb_handle_timing)) call runend('MOD_ALYA2TALP: ERROR WHEN ASSOCIATING HANDLE_TIMING')
    if (.not. c_associated(dlb_handle_global)) call runend('MOD_ALYA2TALP: ERROR WHEN ASSOCIATING HANDLE_GLOBAL')
    do imodu = 0,mmodu
       if (.not. c_associated(dlb_module(imodu) % dlb_handle) ) call runend('MOD_ALYA2TALP: ERROR WHEN ASSOCIATING HANDLE_MODULE')
    end do

    call alya2talp_MonitoringRegionStart(GLOBAL_REGION=.true.,TIMING_REGION=.true.)
#endif
    
  end subroutine alya2talp_register
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Start monitoring
  !> @details Start monitoring a region
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2talp_MonitoringRegionStart(&
       DYNAMIC_ALLOCATION_REGION,&
       GLOBAL_REGION,&
       TIMING_REGION,&
       MODULE_REGION,&
       CURRENT_MODULE)

    logical(lg), optional, intent(in) :: DYNAMIC_ALLOCATION_REGION
    logical(lg), optional, intent(in) :: GLOBAL_REGION
    logical(lg), optional, intent(in) :: TIMING_REGION
    logical(lg), optional, intent(in) :: MODULE_REGION
    integer(ip), optional, intent(in) :: CURRENT_MODULE
    integer(4)                        :: ierr
    logical(lg)                       :: if_dynamic_allocation_region
    logical(lg)                       :: if_timing_region
    logical(lg)                       :: if_global_region
    logical(lg)                       :: if_module_region
    integer(ip)                       :: imodu

#if defined ALYA_TALP
    if( present(DYNAMIC_ALLOCATION_REGION) ) then
       if_dynamic_allocation_region = DYNAMIC_ALLOCATION_REGION
    else
       if_dynamic_allocation_region = .false.
    end if
    if( present(TIMING_REGION) ) then
       if_timing_region = TIMING_REGION
    else
       if_timing_region = .false.
    end if
     if( present(GLOBAL_REGION) ) then
       if_global_region = GLOBAL_REGION
    else
       if_global_region = .false.
    end if
    if( present(MODULE_REGION) ) then
       if_module_region = MODULE_REGION
    else
       if_module_region = .false.
    end if
    if( present(CURRENT_MODULE) ) then
       imodu = CURRENT_MODULE
    else
       imodu = modul
    end if    

    if( if_dynamic_allocation_region ) then
       ierr = DLB_MonitoringRegionStart(dlb_handle_da)
    end if
    if( if_global_region ) then
       ierr = DLB_MonitoringRegionStart(dlb_handle_global)
    end if
    if( if_timing_region ) then
       ierr = DLB_MonitoringRegionStart(dlb_handle_timing)
    end if
    if( if_module_region ) then
       ierr = DLB_MonitoringRegionStart(dlb_module(imodu) % dlb_handle)
    end if
#endif

  end subroutine alya2talp_MonitoringRegionStart

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Stop monitoring
  !> @details Stop monitoring a region
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2talp_MonitoringRegionStop(&
       DYNAMIC_ALLOCATION_REGION,&
       GLOBAL_REGION,&
       TIMING_REGION,&
       MODULE_REGION,&
       CURRENT_MODULE)

    logical(lg), optional, intent(in) :: DYNAMIC_ALLOCATION_REGION
    logical(lg), optional, intent(in) :: GLOBAL_REGION
    logical(lg), optional, intent(in) :: TIMING_REGION
    logical(lg), optional, intent(in) :: MODULE_REGION
    integer(ip), optional, intent(in) :: CURRENT_MODULE
    integer(4)                        :: ierr
    logical(lg)                       :: if_dynamic_allocation_region
    logical(lg)                       :: if_timing_region
    logical(lg)                       :: if_global_region
    logical(lg)                       :: if_module_region
    integer(ip)                       :: imodu

#if defined ALYA_TALP
    if( present(DYNAMIC_ALLOCATION_REGION) ) then
       if_dynamic_allocation_region = DYNAMIC_ALLOCATION_REGION
    else
       if_dynamic_allocation_region = .false.
    end if
    if( present(TIMING_REGION) ) then
       if_timing_region = TIMING_REGION
    else
       if_timing_region = .false.
    end if
     if( present(GLOBAL_REGION) ) then
       if_global_region = GLOBAL_REGION
    else
       if_global_region = .false.
    end if
   if( present(MODULE_REGION) ) then
       if_module_region = MODULE_REGION
    else
       if_module_region = .false.
    end if
    if( present(CURRENT_MODULE) ) then
       imodu =  CURRENT_MODULE
    else
       imodu = modul
    end if    

    if( if_dynamic_allocation_region ) then
       ierr = DLB_MonitoringRegionStop(dlb_handle_da)
    end if
    if( if_global_region ) then
       ierr = DLB_MonitoringRegionStop(dlb_handle_global)
    end if
    if( if_timing_region ) then
       ierr = DLB_MonitoringRegionStop(dlb_handle_timing)
    end if
    if( if_module_region ) then
       ierr = DLB_MonitoringRegionStop(dlb_module(imodu) % dlb_handle)
    end if
#endif

  end subroutine alya2talp_MonitoringRegionStop

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Reset a region
  !> @details Reset a region: put counter to zero
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2talp_MonitoringRegionReset(&
       DYNAMIC_ALLOCATION_REGION,&
       GLOBAL_REGION,&
       TIMING_REGION,&
       MODULE_REGION,&
       CURRENT_MODULE)

    logical(lg), optional, intent(in) :: DYNAMIC_ALLOCATION_REGION
    logical(lg), optional, intent(in) :: GLOBAL_REGION
    logical(lg), optional, intent(in) :: TIMING_REGION
    logical(lg), optional, intent(in) :: MODULE_REGION
    integer(ip), optional, intent(in) :: CURRENT_MODULE
    integer(4)                        :: ierr
    logical(lg)                       :: if_dynamic_allocation_region
    logical(lg)                       :: if_global_region
    logical(lg)                       :: if_timing_region
    logical(lg)                       :: if_module_region
    integer(ip)                       :: imodu

#if defined ALYA_TALP
    if( present(DYNAMIC_ALLOCATION_REGION) ) then
       if_dynamic_allocation_region = DYNAMIC_ALLOCATION_REGION
    else
       if_dynamic_allocation_region = .false.
    end if
    if( present(GLOBAL_REGION) ) then
       if_global_region = GLOBAL_REGION
    else
       if_global_region = .false.
    end if
    if( present(TIMING_REGION) ) then
       if_timing_region = TIMING_REGION
    else
       if_timing_region = .false.
    end if
    if( present(MODULE_REGION) ) then
       if_module_region = MODULE_REGION
    else
       if_module_region = .false.
    end if

    if( if_dynamic_allocation_region ) then
       ierr = DLB_MonitoringRegionreset(dlb_handle_da)
    end if
    if( if_global_region ) then
       ierr = DLB_MonitoringRegionreset(dlb_handle_global)
    end if
    if( if_timing_region ) then
       ierr = DLB_MonitoringRegionreset(dlb_handle_timing)
    end if
    if( if_module_region ) then
       do imodu = 0,mmodu
          if( kfl_modul(imodu) == 1 ) &
               ierr = DLB_MonitoringRegionreset(dlb_module(imodu) % dlb_handle)
       end do
    end if
#endif

  end subroutine alya2talp_MonitoringRegionReset

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Get values from a region
  !> @details Get TALP counter from a region. Time is given in ns
  !>          by TALP.
  !> 
  !-----------------------------------------------------------------------

#if defined ALYA_TALP
  subroutine alya2talp_time(elapsed_time,accumulated_MPI_time,accumulated_computation_time,dlb_handle,dlb_monitor)
    
    real(rp),            optional, intent(out) :: elapsed_time
    real(rp),            optional, intent(out) :: accumulated_MPI_time
    real(rp),            optional, intent(out) :: accumulated_computation_time
    type(c_ptr),         optional              :: dlb_handle
    type(dlb_monitor_t), optional, pointer     :: dlb_monitor

    if( present(dlb_handle) .and. present(dlb_monitor) ) then
       call c_f_pointer(dlb_handle, dlb_monitor)
       if( present(elapsed_time) )                 elapsed_time                 = real(dlb_monitor % elapsed_time,rp)*1.0e-9_rp
       if( present(accumulated_MPI_time) )         accumulated_MPI_time         = real(dlb_monitor % accumulated_MPI_time,rp)*1.0e-9_rp
       if( present(accumulated_computation_time) ) accumulated_computation_time = real(dlb_monitor % accumulated_computation_time,rp)*1.0e-9_rp
    else
       call c_f_pointer(dlb_handle_da, dlb_monitor_da)
       if( present(elapsed_time) )                 elapsed_time                 = real(dlb_monitor_da % elapsed_time,rp)*1.0e-9_rp
       if( present(accumulated_MPI_time) )         accumulated_MPI_time         = real(dlb_monitor_da % accumulated_MPI_time,rp)*1.0e-9_rp
       if( present(accumulated_computation_time) ) accumulated_computation_time = real(dlb_monitor_da % accumulated_computation_time,rp)*1.0e-9_rp
    end if

  end subroutine alya2talp_time
#endif

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Get values from a region
  !> @details Get TALP counter from a region. Time is given in ns
  !>          by TALP.
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2talp_elapsed_time(elapsed_time)
    
    real(rp), intent(out) :: elapsed_time

#if defined ALYA_TALP
    call c_f_pointer(dlb_handle_timing, dlb_monitor_timing)
    elapsed_time = real(dlb_monitor_timing % elapsed_time,rp)*1.0e-9_rp
#endif

  end subroutine alya2talp_elapsed_time

  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2020-03-06
  !> @brief   Performance indicators
  !> @details Parallel peformance indicators
  !>
  !>            Definitions taken from PoP project: https://pop-coe.eu/node/69 
  !>                                                                           
  !>                t_i^w      t_i^MPI                                         
  !>            <------------><------->                                        
  !>                                                                           
  !>            +---------------------+  -+                                    
  !>            |\\\\\\\\\\\\\        +   |                                    
  !>            +---------------------+   |                                    
  !>            |\\\\\\\\\\\\\\\\\    +   | P processors                       
  !>            +---------------------+   |                                    
  !>            |\\\\\\\              +   |                                    
  !>            +---------------------+  -+                                    
  !>                                                                           
  !>            <---------------------> te                                     
  !>                                                                           
  !>            t_i^w ..... working time of process i                          
  !>            t_i^MPI ... MPI communication time of process i                
  !>            t_e ....... Elpased time (same for all)                        
  !>                                                                           
  !>            t_ave^w = sum_i t_i^w / P (average working time)               
  !>            t_max^w = max_i t_i^w     (max working time)                   
  !>                                                                           
  !>                                                                           
  !>            Load balance:             LB = t_ave^w /t_max^w                
  !>            Communication efficiency: CE = t_max^w /t_e                    
  !>            Parallel efficiency:      PE = LB * CE = t_ave^w / te          
  !>
  !-----------------------------------------------------------------------
  
  subroutine alya2talp_parallel_efficiency(&
       PE,&
       LB,&
       CE,&
       time_comp_ave,&
       time_mpi_ave,&
       time_comp_max,&
       time_mpi_max,&
       time_total,&
       GLOBAL_REGION,&
       DYNAMIC_ALLOCATION_REGION,&
       MODULE_REGION,&
       CURRENT_MODULE)
    
    real(rp),              intent(out) :: PE                         !< Parallel efficiency
    real(rp),              intent(out) :: LB                         !< Load balance
    real(rp),              intent(out) :: CE                         !< Communication efficiency
    real(rp),    optional, intent(out) :: time_comp_ave
    real(rp),    optional, intent(out) :: time_mpi_ave
    real(rp),    optional, intent(out) :: time_comp_max
    real(rp),    optional, intent(out) :: time_mpi_max
    real(rp),    optional, intent(out) :: time_total
    logical(lg), optional, intent(in)  :: GLOBAL_REGION
    logical(lg), optional, intent(in)  :: DYNAMIC_ALLOCATION_REGION
    logical(lg), optional, intent(in)  :: MODULE_REGION
    integer(ip), optional, intent(in)  :: CURRENT_MODULE
    integer(ip)                        :: imodu
    logical(lg)                        :: if_global_region
    logical(lg)                        :: if_dynamic_allocation_region
    logical(lg)                        :: if_module_region
    real(rp)                           :: elapsed_time
    real(rp)                           :: accumulated_MPI_time
    real(rp)                           :: accumulated_computation_time
    real(rp)                           :: ave_elapsed_time
    real(rp)                           :: max_elapsed_time
    real(rp)                           :: max_accumulated_MPI_time
    real(rp)                           :: ave_accumulated_MPI_time
    real(rp)                           :: ave_accumulated_computation_time
    real(rp)                           :: max_accumulated_computation_time
    real(rp)                           :: P

    if( present(GLOBAL_REGION) ) then
       if_global_region = GLOBAL_REGION
    else
       if_global_region = .false.
    end if
    if( present(DYNAMIC_ALLOCATION_REGION) ) then
       if_dynamic_allocation_region = DYNAMIC_ALLOCATION_REGION
    else
       if_dynamic_allocation_region = .false.
    end if
    if( present(MODULE_REGION) ) then
       if_module_region = MODULE_REGION
    else
       if_module_region = .false.
    end if
    if( present(CURRENT_MODULE) ) then
       imodu = CURRENT_MODULE
    else
       imodu = modul
    end if    

#if defined ALYA_TALP

    if( if_global_region ) then
       call alya2talp_time(elapsed_time,accumulated_MPI_time,accumulated_computation_time,&
            dlb_handle_global,&
            dlb_monitor_global)
    else if( if_dynamic_allocation_region ) then
       call alya2talp_time(elapsed_time,accumulated_MPI_time,accumulated_computation_time,&
            dlb_handle_da,&
            dlb_monitor_da)       
    else if( if_module_region ) then
       call alya2talp_time(elapsed_time,accumulated_MPI_time,accumulated_computation_time,&
            dlb_module(imodu) % dlb_handle,&
            dlb_module(imodu) % dlb_monitor)
    else
       call runend('MOD_ALYA2TAP: DO NOT KNOW WHAT TO DO')
    end if

     ave_accumulated_MPI_time         = accumulated_MPI_time
     max_accumulated_MPI_time         = accumulated_MPI_time
     ave_accumulated_computation_time = accumulated_computation_time
     max_accumulated_computation_time = accumulated_computation_time
     ave_elapsed_time                 = elapsed_time
     max_elapsed_time                 = elapsed_time

     call PAR_AVERAGE(ave_accumulated_MPI_time)
     call PAR_MAX    (max_accumulated_MPI_time)
     call PAR_AVERAGE(ave_accumulated_computation_time)
     call PAR_MAX    (max_accumulated_computation_time)
     call PAR_AVERAGE(ave_elapsed_time)
     call PAR_MAX    (max_elapsed_time)
     
     LB = ave_accumulated_computation_time / (max_accumulated_computation_time+zeror) ! Load balance
     CE = max_accumulated_computation_time / (max_elapsed_time+zeror)                 ! Communication efficiency
     PE = ave_accumulated_computation_time / (max_elapsed_time+zeror)                 ! Parallel efficiency

     if( present(time_comp_ave)  ) time_comp_ave  = ave_accumulated_computation_time
     if( present(time_comp_max)  ) time_comp_max  = max_accumulated_computation_time
     if( present(time_mpi_ave)   ) time_mpi_ave   = ave_accumulated_MPI_time
     if( present(time_mpi_max)   ) time_mpi_max   = max_accumulated_MPI_time
     if( present(time_total)     ) time_total     = max_elapsed_time

#else

     PE = -1.0_rp
     LB = -1.0_rp
     CE = -1.0_rp
     if( present(time_comp_ave)  ) time_comp_ave  = -1.0_rp
     if( present(time_comp_max)  ) time_comp_max  = -1.0_rp
     if( present(time_mpi_ave)   ) time_mpi_ave   = -1.0_rp
     if( present(time_mpi_max)   ) time_mpi_max   = -1.0_rp
     if( present(time_total)     ) time_total     = -1.0_rp

#endif

  end subroutine alya2talp_parallel_efficiency

end module mod_alya2talp
!> @}
