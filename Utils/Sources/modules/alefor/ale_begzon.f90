subroutine ale_begzon
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_begzon
  ! NAME 
  !    ale_begste
  ! DESCRIPTION
  !    This routine prepares for a new coupling iteration of the ALE formulation
  !    equation      
  ! USES
  !    ale_updunk
  ! USED BY
  !    Alefor
  !***
  !-----------------------------------------------------------------------
  use def_coupli,    only : kfl_gozon
  use def_kintyp,    only : ip
  implicit none
  integer(ip) :: idime
  !
  ! At the beginning of a coupling iteration, the coordinates  
  ! must go back to the previous time step values in FSI
  !
  if ( kfl_gozon == 1_ip ) call ale_updunk(6_ip)

end subroutine ale_begzon
