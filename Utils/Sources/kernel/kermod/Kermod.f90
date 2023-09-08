subroutine Kermod(order)
  !-----------------------------------------------------------------------
  !****f* kermod/Kermod
  ! NAME 
  !    Kermod
  ! DESCRIPTION
  !    This routine deals with the incompressible NS equations.
  !    Kermod is monitorized for Paraver. 
  ! USES
  !    ker_turnon
  !    ker_timste
  !    ker_begste
  !    ker_doiter
  !    ker_concon
  !    ker_conblk
  !    ker_newmsh
  !    ker_endste
  !    ker_turnof
  ! USED BY
  !    Reapro
  !    Turnon
  !    Timste
  !    Begste
  !    Doiter
  !    Concon
  !    Conblk
  !    Newmsh
  !    Endste
  !    Turnof
  !***
  !-----------------------------------------------------------------------
  use def_master
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: order

  modul = mmodu
  call moddef(9_ip)

  select case ( order )

  case( -ITASK_TURNON )
     call livinf(51_ip,' ',modul) 
     call ker_turnon(1_ip)

  case( ITASK_TURNON )
     call ker_turnon(2_ip)

  case( ITASK_TIMSTE ) 
     !call ker_timste()

  case( -ITASK_INIUNK )
     call livinf(53_ip,' ',modul) 
     call ker_iniunk()

  case( ITASK_BEGSTE ) 
     call ker_begste()

  case( ITASK_DOITER )

     !call ker_doiter()
  case(  ITASK_CONCOU )
     !call ker_concou()

  case(  ITASK_CONBLK )
     !call ker_conblk()

  case(  ITASK_NEWMSH )
     !call ker_newmsh()

  case(  ITASK_ENDSTE )
     call ker_endste()

  case( -ITASK_FILTER )
     call ker_filter()

  case( -ITASK_OUTPUT )
     call ker_output()

  case( ITASK_TURNOF )
     !call ker_turnof()

  case( -1_ip )
     !
     ! Compute wall distance
     !
     call ker_walgen(2_ip)
     call ker_walnor(2_ip)

  end select
  !
  ! Coupling
  ! 
  if( order > 1000 ) call ker_plugin(order-1000_ip) 

  modul = 0
  call moddef(9_ip)

end subroutine Kermod
