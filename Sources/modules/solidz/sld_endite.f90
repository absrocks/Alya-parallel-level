!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_endite.f90
!> @author  Mariano Vazquez
!> @date    August, 2006
!>          - Subroutine written
!> @brief   Convergence checks and unknowns update
!> @details
!>          \verbatim
!>          ITASK = 1 The end of an internal iteration
!>                = 2 The end of an external iteration
!>          \endverbatim
!> @todo
!>           <GGU> Revision sld_coupli (obsolete for nastin and solidz)
!>           <GGU> Outputs for reaction forces and displacement should
!>                 be written for (i) converged iterations (ii) different
!>                 output file (not cvgunk)
!> @}
!------------------------------------------------------------------------

subroutine sld_endite(itask)

  use def_kintyp,            only : ip, rp
  use def_master,            only : itinn, modul, intost
  use def_master,            only : ITASK_ENDINN, ITASK_ENDITE
  use def_master,            only : solve_sol
  use def_master,            only : ittim
  use def_domain,            only : kfl_elcoh
  use mod_messages,          only : livinf
  use def_solidz,            only : miinn_sld
  use def_solidz,            only : kfl_volca_sld, kfl_cycle_sld
  use def_solidz,            only : kfl_xfeme_sld, last_iters_sld
  use def_solidz,            only : kfl_cycle_sld
  use mod_sld_energy,        only : sld_updene
  use mod_sld_fe2,           only : fe2_update_vars, fe2_write_profiling
  use mod_sld_cardiac_cycle, only : sld_cardiac_cycle_compute_cavity_volume, sld_manage_cycles

  implicit none

  integer(ip), intent(in) :: itask   !< 1, inner iterations; 2, outer iterations
  integer(ip)             :: maxiter
  character(300)          :: messa

  !
  ! Calculate the volume inside the cavities and windkessel pressure if needed
  !
  if (kfl_volca_sld==1_ip)  then
     call sld_cardiac_cycle_compute_cavity_volume(itask)
     if (kfl_cycle_sld == 1) call sld_manage_cycles(itask)
  end if

  select case(itask)

  case( ITASK_ENDINN )

     !-------------------------------------------------------------------
     !
     !  Compute convergence residual of the internal iteration
     !
     !-------------------------------------------------------------------

     call sld_cvgunk(ITASK_ENDINN)   ! Residual
     call sld_updunk(ITASK_ENDINN)   ! Update: u(n,i,j-1) <-- u(n,i,j)
     !
     ! Live info iterations
     !
     maxiter = solve_sol(1) % miter
     messa = &
          ' (SUBIT: '//trim(intost(itinn(modul)))//'/'//trim(intost(miinn_sld))//' IT: '//trim(intost(last_iters_sld))//'/'//trim(intost(maxiter))//')'
     call livinf(-3_ip,messa,1_ip)
     call livinf(56_ip,' ',modul)
     !
     ! Reaction forces at direchlet nodes printed in the cvgunk
     ! <GGU> It should not be here.
     call sld_residual_force()
     !
     ! Convergence and Timings
     !
     call sld_cvgunk(0_ip)
     !
     ! Matrix output
     !
     call sld_outmat()

  case( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     !  Compute convergence residual of the external iteration
     !
     !-------------------------------------------------------------------

     call livinf(16_ip,' ',itinn(modul))

     call sld_cvgunk(ITASK_ENDITE)  ! Residual || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||)
     call sld_updunk(ITASK_ENDITE)  ! Update unknowns: u(n,i-1,*) <-- u(n,i,*)
     call sld_updene(ITASK_ENDITE)  ! Update energies
     !
     ! Coupling
     !
     call sld_coupli(ITASK_ENDITE)
     !
     ! Cohesive elements and X-FEM
     !
     if (kfl_elcoh > 0 .or. kfl_xfeme_sld > 0) call sld_updcoh(2_ip)
     !
     ! Write MicroPP profiling
     !
     call fe2_write_profiling(ittim)
     !
     ! Update MicroPP internal variables
     !
     call fe2_update_vars()

  end select

end subroutine sld_endite
