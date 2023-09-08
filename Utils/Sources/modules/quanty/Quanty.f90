subroutine Quanty(order)
  !-----------------------------------------------------------------------
  !****f* Quanty/quanty
  ! NAME
  !   Quanty
  ! DESCRIPTION
  !   TDDFT 
  ! USES
  !    qua_turnon
  !    qua_timste
  !    qua_begste
  !    qua_doiter
  !    qua_concon
  !    qua_conblk
  !    qua_newmsh
  !    qua_endste
  !    qua_turnof
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
  implicit none
  integer(ip), intent(in) :: order

  select case (order)

  case(ITASK_TURNON)
     call qua_turnon()
  case(ITASK_INIUNK)
     call qua_iniunk()
  case(ITASK_TIMSTE) 
     call qua_timste()
  case(ITASK_BEGSTE) 
     call qua_begste()
  case(ITASK_DOITER)
     call qua_doiter()
  case(ITASK_CONCOU)
     call qua_concou()
  case(ITASK_CONBLK)
     call qua_conblk()
!  case(ITASK_NEWMSH)
     !if(kfl_modul(modul)/=0) call qua_newmsh()
  case(ITASK_ENDSTE)
     call qua_endste()
  case(ITASK_TURNOF)
     call qua_turnof()

  end select

end subroutine Quanty
