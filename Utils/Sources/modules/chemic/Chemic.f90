subroutine Chemic(order)
  !-----------------------------------------------------------------------
  !****f* chemic/Chemic
  ! NAME
  !   Chemic
  ! DESCRIPTION
  !   This routine deals with the ADS equation. The task done
  !   corresponds to the order given by the master.
  ! USES
  !    chm_turnon
  !    chm_timste
  !    chm_begste
  !    chm_doiter
  !    chm_concon
  !    chm_conblk
  !    chm_newmsh
  !    chm_endste
  !    chm_turnof
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
     call chm_turnon()
  case(ITASK_INIUNK) 
     call chm_iniunk()
   case(ITASK_TIMSTE) 
     call chm_timste()
  case(ITASK_BEGSTE) 
     call chm_begste()
  case(ITASK_DOITER)
     call chm_doiter()
  case(ITASK_CONCOU)
     call chm_concou()
  case(ITASK_CONBLK)
     call chm_conblk()
  case(ITASK_NEWMSH)
     call chm_newmsh()
  case(ITASK_ENDSTE)
     call chm_endste()
  case(ITASK_OUTPUT)
     call chm_output()
  case(ITASK_TURNOF)
     call chm_turnof()

  end select

end subroutine Chemic
