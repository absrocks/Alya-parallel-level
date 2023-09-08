subroutine Immbou(order)
  !-----------------------------------------------------------------------
  !****f* Immbou/immbou
  ! NAME
  !   Immbou
  ! DESCRIPTION
  !   This routine deals with the immbouature equation. The task done
  !   corresponds to the order given by the master.
  ! USES
  !    ibm_turnon
  !    ibm_timste
  !    ibm_begste
  !    ibm_doiter
  !    ibm_concon
  !    ibm_conblk
  !    ibm_newmsh
  !    ibm_endste
  !    ibm_turnof
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
  use def_domain
  implicit none
  integer(ip), intent(in) :: order

  select case (order)

  case(ITASK_TURNON)
     !if( nimbo > 0 ) return
     call ibm_turnon()
  case(ITASK_TIMSTE)
     call ibm_timste()
  case(ITASK_INIUNK)
     call ibm_iniunk()
     !call ibm_blende(1_ip)
  case(ITASK_BEGSTE)
     call ibm_begste()
  case(ITASK_DOITER)
     call ibm_doiter()
     !call ibm_blende(2_ip)
  case(ITASK_CONCOU)
     !call ibm_concou()
  case(ITASK_CONBLK)
     !call ibm_conblk()
  case(ITASK_NEWMSH)
     !call ibm_newmsh()
  case(ITASK_ENDSTE)
     call ibm_endste()
  case(ITASK_OUTPUT)
     call ibm_output()
  case(ITASK_TURNOF)
     !call ibm_turnof()
  end select


end subroutine Immbou
