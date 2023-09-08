!------------------------------------------------------------------------
!> @defgroup Nastin
!> Incompressible - Boussinesq - Low Mach Navier-Stokes equations
!> @{
!> @file    Nastin.f90
!> @date    10/10/1972
!> @author  Guillaume Houzeaux
!> @brief   Incompressible NSI main
!> @details Nastin: incompressible - Low Mach Navier-Stokes equations. Main subroutine
!> @}
!------------------------------------------------------------------------
subroutine Nastin(order)

  use      def_master
  use      def_nastin

  use mod_commdom_nsi,  only: commdom_nsi_plugin

  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: icoup

  select case ( order )

  case( ITASK_TURNON )
     call nsi_turnon()
  case( ITASK_TIMSTE ) 
     call nsi_timste()
  case( ITASK_INIUNK ) 
     call nsi_iniunk()
  case( ITASK_BEGSTE )
     call nsi_begste()
  case( ITASK_DOITER )
     call nsi_doiter()
  case( ITASK_CONCOU )
     call nsi_concou()
  case( ITASK_CONBLK )
     call nsi_conblk()
  case( ITASK_NEWMSH )
     call nsi_newmsh()
  case( ITASK_ENDSTE )
     call nsi_endste()
  case( ITASK_FILTER )
     call nsi_filter()
  case( ITASK_OUTPUT )
     call nsi_output()
  case( ITASK_TURNOF )
     call nsi_turnof()
 case( ITASK_DOOPTI )
!!     call nsi_doopti()
  end select
  !
  ! Coupling
  ! 
  !if(  order > 1000_ip  ) call nsi_fsiexch( 1_ip,order-1000_ip ) ! Compute and send 
  if( order > 1000 ) call nsi_plugin(order-1000_ip) 

  call commdom_nsi_plugin()  !< 2016MAR29 
  
end subroutine Nastin

