subroutine qua_parall(itask)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_parall
  ! NAME
  !    qua_parall
  ! DESCRIPTION
  !    This routine is a bridge to Parall service  
  ! USED BY
  !    qua_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_quanty
  implicit none
  integer(ip), intent(in) :: itask

  if(kfl_paral>=0) then

     select case(itask)

     case(1)
        !
        ! Exchange data read in qua_reaphy, qua_reanut and qua_reaous
        ! always using MPI, even if this is a partition restart
        !
        call qua_sendat(1_ip)
        call qua_sendat(3_ip)

     case(2)
        !
        ! Exchange data read in qua_reabcs
        !
        call qua_sendat(2_ip)

     case(3) 


     case(4)

     end select

  end if

end subroutine qua_parall
