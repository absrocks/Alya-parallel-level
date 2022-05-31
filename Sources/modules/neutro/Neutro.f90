!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    Neutro.f90
!> @date    29/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Generic module
!> @details Generic module
!> @}
!------------------------------------------------------------------------

subroutine Neutro(order)

  use def_master
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: icoup

  select case ( order )

  case( ITASK_TURNON )
     call neu_turnon()
  case( ITASK_TIMSTE ) 
     call neu_timste()
  case( ITASK_INIUNK ) 
     call neu_iniunk()
  case( ITASK_BEGSTE )
     call neu_begste()
  case( ITASK_DOITER )
     call neu_doiter()
  case( ITASK_CONCOU )
     call neu_concou()
  case( ITASK_CONBLK )
     call neu_conblk()
  case( ITASK_NEWMSH )
     ! call neu_newmsh()
  case( ITASK_ENDSTE )
     call neu_endste()
  case( ITASK_FILTER )
     ! call neu_filter()
  case( ITASK_OUTPUT )
     call neu_output()
  case( ITASK_TURNOF )
     call neu_turnof()
  end select
  !
  ! Coupling
  ! 
  !if( order > 1000 ) call neu_plugin(order-1000_ip) 
  
end subroutine Neutro

