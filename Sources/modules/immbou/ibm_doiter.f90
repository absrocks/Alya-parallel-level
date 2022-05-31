subroutine ibm_doiter()
  !-----------------------------------------------------------------------
  !****f* ibm_doiter/ibm_doiter
  ! NAME
  !    ibm_doiter
  ! DESCRIPTION
  !    This routines solves the Euler and Newton equations for rigid bodies
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_immbou
  use mod_messages, only : livinf
  implicit none
  !
  ! Steady state
  !
  if( kfl_timei_ibm == 0 ) return
  !
  ! Begin iteration
  !
  itinn(modul)  = 0
  call livinf( 56_ip,' ',modul)
  call livinf(160_ip,' ',modul)
  !
  ! Solve
  !
  call ibm_solite(1_ip)
  !
  ! Coupling
  !
  ! Call before in ibm_solite
  !
  !call ibm_coupli(ITASK_ENDITE)

end subroutine ibm_doiter
