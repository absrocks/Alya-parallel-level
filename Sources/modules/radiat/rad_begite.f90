
subroutine rad_begite
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_begite
  ! NAME 
  !    rad_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the radiation
  !    equation
  ! USES
  !    rad_tittim
  !    rad_updbcs
  !    rad_inisol
  !    rad_updunk
  ! USED BY
  !    rad_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_radiat
  use mod_messages, only : livinf
  implicit none
  !
  ! Initializations
  !
  kfl_goite_rad = 1 
  itinn(modul)  = 0
  if(itcou==1) call rad_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Update boundary conditions
  !
  call rad_updbcs(two)
  !
  ! Set up the solver parameters for the radiation equation
  !
  call rad_inisol()
  !
  ! Obtain the initial guess for inner iterations
  !
  call rad_updunk(two)

end subroutine rad_begite
    
