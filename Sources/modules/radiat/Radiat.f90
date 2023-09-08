subroutine Radiat(order)
  !-----------------------------------------------------------------------
  !****f* Radiat
  ! NAME
  !   Radiat
  ! DESCRIPTION
  !   This module deals with radiation heat transfer. The task done
  !   corresponds to the order given by the master.
  ! USES
  !    rad_turnon
  !    rad_timste
  !    rad_begste
  !    rad_doiter
  !    rad_concon
  !    rad_conblk
  !    rad_newmsh
  !    rad_endste
  !    rad_turnof
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
  use def_solver
  implicit none
  integer(ip), intent(in) :: order

  select case (order) 

  case(ITASK_TURNON)
     call rad_turnon()  !!!F Details left (maybe works)
  case(ITASK_TIMSTE) 
     call rad_timste()  !!F Someday time dependence will enter this
  case(ITASK_INIUNK) 
     call rad_iniunk()  
  case(ITASK_BEGSTE) 
     call rad_begste()  !!F Need to update BCS and check from previous step(maybe)
  case(ITASK_DOITER)
     call rad_doiter()  
  case(ITASK_CONCOU)
     call rad_concou()  
  case(ITASK_CONBLK)
     call rad_conblk() 
  case(ITASK_NEWMSH)
     call rad_newmsh() 
  case(ITASK_ENDSTE)
     call rad_endste() 
  case(ITASK_OUTPUT)
     call rad_output()  !!F Check how sets and witnesses work...
  case(ITASK_TURNOF)
     call rad_turnof()  ! Redo Latex output

  end select

end subroutine Radiat
