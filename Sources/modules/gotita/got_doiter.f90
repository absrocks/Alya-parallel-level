subroutine got_doiter
!-----------------------------------------------------------------------
!****f* Gotita/got_doiter
! NAME 
!    got_doiter
! DESCRIPTION
!    This routine solves an iteration of the linearized incomcdropible NS
!    equations.
! USES
!    got_begite
!    got_updbcs
!    got_solite
! USED BY
!    Gotita
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_solver
  use      def_gotita
  implicit none
   
  if(kfl_stead_got==0) then
     call got_begite()
     do while(kfl_goite_got==1)
        call got_solite()
        call got_endite(one)
     end do
     call got_endite(two)
  end if

end subroutine got_doiter

