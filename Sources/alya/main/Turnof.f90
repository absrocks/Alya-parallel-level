!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    Turnof.f90
!> @author  houzeaux
!> @date    2019-02-27
!> @brief   Turnof 
!> @details Finish Alya... smoothly...
!> @} 
!-----------------------------------------------------------------------

subroutine Turnof
  
  use def_kintyp,               only : ip
  use def_master,               only : IPARALL
  use def_master,               only : ITASK_TURNOF
  use mod_par_output_partition, only : par_output_partition
  use mod_par_output_partition, only : par_output_solvers
  use mod_cou_output,           only : cou_output_timings
  use mod_bourgogne_pinotnoir,  only : bourgogne
  use mod_mpio_par_configure,   only : par_mpio_finalize
  use mod_messages,             only : messages_live
  use mod_moduls,               only : moduls
  use mod_repartitioning,       only : repartitioning_timings
  use mod_performance,          only : performance_outcpu
  use mod_performance,          only : performance_outmem
  use mod_performance,          only : performance_output
  use mod_parall_destructor,    only : parall_destructor
  use mod_mpio_par_postpr,      only : posmpio_destructor
  use mod_mass_matrix,          only : mass_matrix_destructor
  use mod_materials,            only : materials_destructor
  use mod_exterior_normal,      only : exterior_normal_destructor
  implicit none
  !
  ! NINJAAAAAAAA
  !
#ifdef NINJA
  call lasterrorninja()
#endif
  call bourgogne(1_ip)
  !
  ! Write info about direct solvers
  !
  call out_direct_solvers()
  !
  ! Writes memory used
  !
  call performance_outmem()
  !
  ! Repartitioning
  !
  call repartitioning_timings()     
  !
  ! Write info about timings in partition file
  !
  call par_output_partition(2_ip)
  !
  ! Postprocess solver information
  !
  call par_output_solvers()
  !
  ! Write latex info file
  !
  call outlat(2_ip)
  !
  ! Turn off modules
  !
  call moduls(ITASK_TURNOF)
  !
  ! Coupling CPU
  !
  call cou_output_timings()
  !
  ! Write info about modules and close module files
  !
  call moddef(4_ip)
  !
  ! Write CPU time heading and master's CPU time
  !
  call performance_outcpu()
  !
  ! Write performance file in csv format
  !
  call performance_output()
  !
  ! End-up latex file
  !
  call outlat(3_ip)
  if (IPARALL) then
    call par_mpio_finalize()
  end if
  call par_turnof()    
  !
  ! moving vtk data
  !
  !call vtkmov(1_ip)
  !
  !call finalize coprocessing
  !
#ifdef CATA
  call messages_live('CALL FINALIZE COPROCESSING')
  call coprocessorfinalize()
#endif

  call bourgogne(2_ip)
  !
  ! Destroy
  !
  call coupli_destructor()
  call cou_memory(0_ip)
  call cshder(2_ip)
  call solver_destructor() 
  call domain_destructor()
  call geometry_destructor()
  call mass_matrix_destructor()
  call exterior_normal_destructor()
  call materials_destructor()
  call parall_destructor()
  call posmpio_destructor()
  !
  ! Uncomment to check state of memory
  !
  !call performance_outmem()
  !
  ! Stop the run
  !
  call runend('O.K.!')
  
end subroutine Turnof
