subroutine Wavequ(order)
  !-----------------------------------------------------------------------
  !****f* Wavequ/wavequ
  ! NAME
  !   Wavequ
  ! DESCRIPTION
  !   This routine deals with the wave equation. The task done
  !   corresponds to the order given by the master.
  !
  !   Created:     20 January 2007 
  !   by:          Guillaume Houzeaux
  !   Responsible: Anne-Cecile Lesage  
  ! USES
  !    wav_turnon
  !    wav_timste
  !    wav_begste
  !    wav_doiter
  !    wav_concon
  !    wav_conblk
  !    wav_newmsh
  !    wav_endste
  !    wav_turnof
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
  use      def_master
  implicit none
  integer(ip), intent(in) :: order

  select case (order)

  case(ITASK_TURNON)
     call wav_turnon()
  case(ITASK_TIMSTE) 
     call wav_timste()
  case(ITASK_INIUNK) 
     call wav_iniunk()
  case(ITASK_BEGSTE) 
     call wav_begste()
  case(ITASK_DOITER)
     call wav_doiter()
  case(ITASK_CONCOU)
     call wav_concou()
  case(ITASK_CONBLK)
     call wav_conblk()
  case(ITASK_NEWMSH)
     call wav_newmsh()
  case(ITASK_ENDSTE)
     call wav_endste()
  case(ITASK_OUTPUT)
     call wav_output()
  case(ITASK_TURNOF)
     call wav_turnof()

  end select


end subroutine Wavequ
