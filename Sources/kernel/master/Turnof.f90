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
  use mod_messages,             only : livinf
  implicit none

  call bourgogne(1_ip)
  !
  ! Live information
  !
  call livinf(10_ip,' ',0_ip)
  !
  ! Write info about direct solvers
  !
  call out_direct_solvers()
  !
  ! Writes memory used
  !
  call outmem()
  !
  ! Write CPU time heading and master's CPU time
  !
  call outcpu()
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
  ! End-up latex file
  !
  call outlat(3_ip)
  !
  ! Turn off Services
  !
  call Hdfpos(10_ip)        ! HDFPOS
  if (IPARALL) then
    call par_mpio_finalize()
  end if
  call Parall( 7_ip)        ! PARALL
  call Dodeme(12_ip)        ! DODEME
  !
  ! moving vtk data
  !
  !call vtkmov(1_ip)
  !
  ! Close postprocess files
  !
  call openfi(4_ip)
  call openfi(6_ip)
  !
  ! Deallocate memory
  !
  call deaker()
  !
  !call finalize coprocessing
  !
#ifdef CATA
  call livinf(204_ip,' ',0_ip)
  call coprocessorfinalize()
#endif

  call bourgogne(2_ip)
  !
  ! Stop the run
  !
  call runend('O.K.!')
  !
end subroutine Turnof
