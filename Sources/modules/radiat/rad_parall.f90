subroutine rad_parall(itask)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_parall
  ! NAME
  !    rad_parall
  ! DESCRIPTION
  !    This routine is a bridge to Parall service  
  ! USED BY
  !    rad_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_radiat
  implicit none
  integer(ip), intent(in) :: itask

  if(kfl_paral>=0) then

     select case(itask)

     case(1)
        !
        ! Exchange data read in rad_reaphy, rad_reanut and rad_reaous
        ! always using MPI, even if this is a partition restart
        !
        call rad_sendat(1_ip)

     case(2)
        !
        ! Exchange data read in rad_reabcs
        !
        call rad_sendat(2_ip)

     case(3) 


     case(4)

     end select

  end if

end subroutine rad_parall
