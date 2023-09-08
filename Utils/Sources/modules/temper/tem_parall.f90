subroutine tem_parall(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_parall
  ! NAME
  !    tem_parall
  ! DESCRIPTION
  !    This routine is a bridge to Parall service  
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_temper
  implicit none
  integer(ip), intent(in) :: itask

  if(kfl_paral>=0) then

     select case(itask)

     case(1)
        !
        ! Exchange data read in tem_reaphy, tem_reanut and tem_reaous
        ! always using MPI, even if this is a partition restart
        !
        call tem_sendat(1_ip)

     case(2)
        !
        ! Exchange data read in tem_reabcs
        !
        call tem_sendat(2_ip)

     case(3) 


     case(4)

     end select

  end if

end subroutine tem_parall
