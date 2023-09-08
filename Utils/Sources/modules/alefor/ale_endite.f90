subroutine ale_endite()
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_endite
  ! NAME 
  !    ale_endite
  ! DESCRIPTION
  !    End of iterations   
  ! USES
  !    ale_coupli
  ! USED BY
  !    ale_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_alefor
  implicit none
  !
  ! Coupling
  !
  call ale_coupli(ITASK_ENDITE)
  !
  ! Update:   VELOL(:,2) = VELOL(:,1)
  !
  call ale_updunk(24_ip)

end subroutine ale_endite

