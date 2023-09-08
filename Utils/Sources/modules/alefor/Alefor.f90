subroutine Alefor(order)
!-----------------------------------------------------------------------
!****f* Alefor/alefor
! NAME
!   Alefor
! DESCRIPTION
!   This routine deals with the ALE formulation equation. The task done
!   corresponds to the order given by the master.
! USES
!    ale_turnon
!    ale_begste
!    ale_doiter
!    ale_gencon
!    ale_endste
!    ale_turnof
! USED BY
!    Turnon
!    Begste
!    Doiter
!    Gencon
!    Endste
!    Turnof
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master

  use mod_commdom_ale,  only: commdom_ale_plugin

  implicit none
  integer(ip) :: order

  select case (order)

  case(ITASK_TURNON)
     call ale_turnon()
  case(ITASK_TIMSTE) 
     !call ale_timste()
  case(ITASK_INIUNK) 
     call ale_iniunk()
  case(ITASK_BEGSTE) 
     call ale_begste()
  case(ITASK_BEGZON)
     call ale_begzon()
  case(ITASK_DOITER)
     call ale_doiter()
  case(ITASK_CONCOU)
     call ale_concou()
  case(ITASK_CONBLK)
     !call ale_conblk()
  case(ITASK_NEWMSH)
     !call ale_newmsh()
  case(ITASK_ENDSTE)
     call ale_endste()
  case(ITASK_OUTPUT)
     call ale_output()
  case(ITASK_TURNOF)
     call ale_turnof()
  end select
  !
  ! Coupling
  !
  if( order > 1000_ip ) call ale_plugin(order-1000) ! Compute and send 

  call commdom_ale_plugin() !< 2016MAR29  
 
end subroutine Alefor
      
