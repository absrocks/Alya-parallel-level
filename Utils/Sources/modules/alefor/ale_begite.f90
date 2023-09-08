subroutine ale_begite
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_begite
  ! NAME 
  !    ale_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the ALE formulation
  !    equation      
  ! USES
  !    ale_inisol
  !    ale_updunk
  ! USED BY
  !    tem_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_alefor
  implicit none
  !
  ! Update boundary conditions
  !
  call ale_updbcs()
  !
  ! Obtain the initial guess for inner iterations
  !
  call ale_updunk(2_ip)

end subroutine ale_begite

