subroutine got_begite
  !-----------------------------------------------------------------------
  !****f* Gotita/got_begite
  ! NAME 
  !    got_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the incomcdropible NS
  !    equations. 
  ! USES
  !    got_tittim
  !    sni_updbcs
  !    got_inisol
  !    got_updunk
  ! USED BY
  !    got_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_gotita
  use mod_messages, only : livinf
  implicit none
  !
  ! Initializations
  !
  kfl_goite_got = 1
  itinn(modul)  = 0
  if(miinn_got==0) kfl_goite_got=0
  if(itcou==1)     call got_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Set up the solver parameters for the NS equations
  !
  call got_inisol(1_ip)
  !
  ! Obtain the initial guess for inner iterations, with the prescriptions
  ! in local coordinates.
  !
  call got_updunk(2_ip)

end subroutine got_begite
