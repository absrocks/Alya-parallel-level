 !-----------------------------------------------------------------------
!> @addtogroup Initia
!> @{
!> @file    Turnon.f90
!> @author  Guillaume Houzeaux
!> @brief   Turnon the run
!> @details Initial operatons: read general data, mesh and modules data
!>          Construct the mesh dependent arrays.
!>
!> @}
!-----------------------------------------------------------------------
subroutine Initia

  use def_parame
  use def_elmtyp
  use def_master
  use def_inpout
  use def_domain
  use def_kermod
  use mod_bourgogne_pinotnoir
  use mod_unity_tests
  use mod_finite_volume,             only : finite_volume_arrays
  use mod_auto_tuning,               only : auto_tuning_SpMV_OpenMP
  use mod_messages,                  only : messages_report
  use mod_communications,            only : PAR_BARRIER
  use mod_par_affinity,              only : par_affinity
  use mod_messages,                  only : messages_live
  use mod_outfor,                    only : outfor
  use mod_messages,                  only : livinf
  use mod_parall_openmp,             only : parall_openmp_adjacency_ompss_unity_test
  use mod_alya2dlb,                  only : alya2dlb_initialization
  use mod_alya2talp,                 only : alya2talp_register
  use mod_alya2signal,               only : alya2signal
  use mod_optimum_partition,         only : optimum_partition_setup
  implicit none

  integer(8) :: count_rate8
 
  !---------------------------------------------------------------------------------
  !
  ! To force Nanos to work with any MPI versions... here just in case!
  !
  !---------------------------------------------------------------------------------

#ifdef ALYA_OMPSS
  integer :: i
  !$omp do schedule(static)
  do i=1, 2
  end do
  !$omp do schedule(dynamic)
  do i=1, 2
  end do
#endif  

  !---------------------------------------------------------------------------------
  !
  ! Start run, creates code communicators in case of coupling, and read problem data
  !
  !---------------------------------------------------------------------------------

  call bourgogne(1_ip)
  !
  ! Compute time rate
  !
  call system_clock(count_rate=count_rate8)
  rate_time = 1.0_rp / max(real(count_rate8,rp),zeror)
  !
  ! Splits the MPI_COMM_WORLD for coupling with other codes. This defines the Alya world
  !
  call par_code_split_universe()
  !
  ! Splits the MPI_COMM_WORLD for coupling
  !
  call par_code_split_world()
  !
  ! Initialize DLB, just after initializing MPI
  !
  call alya2dlb_initialization()  
  call alya2talp_register()
  !
  ! Initializations
  !
  call inirun()
  !
  ! Read problem data
  !
  call Reapro()
  !
  ! Initialize time counters
  !
  call setgts(ITASK_INITIA)
  !
  ! Initializations
  !
  call inidom()
  !
  ! Define element types
  !
  call elmtyp()                            ! Element lists: NNODE, LTOPO, LDIME, LLAPL...
  !
  ! Open module data files
  !
  call moddef(1_ip)
  !
  ! Affinity of processes
  !
  call par_affinity()
  !
  ! Signal handling
  !
  call alya2signal()
  !
  ! Automatic partitioning with Pycomms
  !
  call optimum_partition_setup()

end subroutine Initia
