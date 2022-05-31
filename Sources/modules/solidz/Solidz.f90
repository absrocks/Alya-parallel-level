!------------------------------------------------------------------------
!> @defgroup Solidz
!> @{
!> @file    Solidz.f90
!> @date    11/07/2018
!> @author  Guillaume Houzeaux
!> @brief   Solidz : Solid mechanics equations. Main subroutine
!> @details Solidz : Solid mechanics equations. Main subroutine
!> @}
!------------------------------------------------------------------------

subroutine Solidz(order)

  use def_kintyp,          only : ip
  use def_master,          only : ITASK_TURNON, ITASK_TURNOF
  use def_master,          only : ITASK_BEGSTE, ITASK_ENDSTE
  use def_master,          only : ITASK_INIUNK, ITASK_DOITER
  use def_master,          only : ITASK_TIMSTE, ITASK_OUTPUT
  use def_master,          only : ITASK_CONCOU, ITASK_CONBLK
  use def_master,          only : ITASK_SOLMEM, ITASK_REDIST
  use def_master,          only : ITASK_BEGRUN
  use def_master,          only : ITASK_READ_RESTART, ITASK_WRITE_RESTART
#ifdef COMMDOM
  use mod_sld_commdom,     only : commdom_sld_plugin
#endif

  implicit none

  integer(ip), intent(in) :: order

  select case (order)

  case( ITASK_TURNON )
     call sld_turnon()
  case( ITASK_BEGRUN )
     call sld_begrun()
  case( ITASK_SOLMEM )
     call sld_solmem()
  case( ITASK_TIMSTE )
     call sld_timste()
  case( ITASK_INIUNK )
     call sld_iniunk()
  case( ITASK_BEGSTE )
     call sld_begste()
  case( ITASK_DOITER )
     call sld_doiter()
  case( ITASK_CONCOU )
     call sld_concou()
  case( ITASK_CONBLK )
     call sld_conblk()
  case( ITASK_ENDSTE )
     call sld_endste()
  case( ITASK_OUTPUT )
     call sld_output()
  case( ITASK_TURNOF )
     call sld_turnof()
  case( ITASK_REDIST )
     call sld_redist()
  case( ITASK_READ_RESTART )
     call sld_restar(ITASK_READ_RESTART)
  case( ITASK_WRITE_RESTART )
     call sld_restar(ITASK_WRITE_RESTART)

  end select
  !
  ! Couplings
  !
  if ( order > 1000_ip ) call sld_plugin(order-1000_ip) ! Compute and send
#ifdef COMMDOM
  call commdom_sld_plugin()
#endif

end subroutine Solidz
