subroutine Gotita(order)
  !-----------------------------------------------------------------------
  !****f* gotita/Gotita
  ! NAME 
  !    Gotita
  ! DESCRIPTION
  !    This routine deals with the incomcdropible NS equations.
  !    Gotita is monitorized for Paraver. The events and counters:
  !    - got_elmope: 100
  !       - got_elmgat: 102
  !       - got_elmtss: 104
  !       - got_elmpre: 106
  !       - got_elmpro: 108
  !       - got_elmres: 110
  !       - got_elmrhs: 112
  !       - got_elmsgs: 114
  !       - got_elmtes: 116
  !       - got_elmmom: 118
  !       - got_elmcon: 120
  !       - got_elmmat: 122
  !       - got_elmdir: 124 
  !       - got_assres: 126
  !       - got_assrhs: 128
  !       - got_assmat: 130
  !    - got_solite: 200
  ! USES
  !    got_turnon
  !    got_timste
  !    got_begste
  !    got_doiter
  !    got_concon
  !    got_conblk
  !    got_newmsh
  !    got_endste
  !    got_turnof
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
     call got_turnon()
  case(ITASK_INIUNK)
     call got_iniunk()
  case(ITASK_TIMSTE) 
     call got_timste()
  case(ITASK_BEGSTE) 
     call got_begste()
  case(ITASK_DOITER)
     call got_doiter()
  case(ITASK_CONCOU)
     call got_concou()
  case(ITASK_CONBLK)
     call got_conblk()
  case(ITASK_NEWMSH)
     call got_newmsh()
  case(ITASK_ENDSTE)
     call got_endste()
  case(ITASK_OUTPUT)
     call got_output()
  case(ITASK_TURNOF)
     call got_turnof()

  end select

end subroutine Gotita
