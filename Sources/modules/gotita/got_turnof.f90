subroutine got_turnof
  !------------------------------------------------------------------------
  !****f* Gotita/got_turnof
  ! NAME 
  !    got_turnof
  ! DESCRIPTION
  !    This routine closes GOTITA module
  ! USES
  !    got_output
  !    got_outlat
  !    got_openfi
  ! USED BY
  !    Gotita
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_gotita
  use      def_solver
  implicit none
  !
  ! Output latex file
  !
  call got_outlat(two)

end subroutine got_turnof

