subroutine got_parall(itask)
!-----------------------------------------------------------------------
!****f* Gotita/got_parall
! NAME
!    got_parall
! DESCRIPTION
!    This routine is a bridge to Parall service  
! USED BY
!    Gotita
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_gotita
  use      mod_memchk
  implicit none
  integer(ip), intent(in) :: itask

  if (kfl_paral>=0) then

     select case (itask)
        
     case(1)
        call got_sendat(one)
     case(2)
        call got_sendat(two)
     case(3)
        
     end select
     
  end if

end subroutine got_parall
