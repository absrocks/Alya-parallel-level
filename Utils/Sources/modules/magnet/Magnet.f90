subroutine Magnet(order)
!-----------------------------------------------------------------------
!****f* Magnet/Magnet
! NAME 
!    Magnet
! DESCRIPTION
!    This routine deals with a magneto dynamics problem 
!    (no displacement currents are present)
! USES
!    mag_turnon
!    mag_begste
!    mag_doiter
!    mag_gencon
!    mag_endste
!    mag_turnof
! USED BY
!    Turnon
!    Begste
!    Doiter
!    Gencon
!    Endste
!    Turnof
!***
!
!    The task done corresponds to the order given by the master.
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  implicit none
  integer(ip), intent(in) :: order

  select case (order)

  case(ITASK_TURNON)
!     if(kfl_modul(modul)/=0) call mag_turnon
  case(ITASK_TIMSTE) 
!     if(kfl_modul(modul)/=0) call mag_timste
  case(ITASK_BEGSTE) 
!     if(kfl_modul(modul)/=0) call mag_begste
  case(ITASK_DOITER)
!     if(kfl_modul(modul)/=0) call mag_doiter
  case(ITASK_CONCOU)
!     if(kfl_modul(modul)/=0) call mag_concou
  case(ITASK_CONBLK)
!     if(kfl_modul(modul)/=0) call mag_conblk
  case(ITASK_NEWMSH)
!!!!     if(kfl_modul(modul)/=0) call mag_newmsh
  case(ITASK_ENDSTE)
!     if(kfl_modul(modul)/=0) call mag_endste
  case(ITASK_TURNOF)
!     if(kfl_modul(modul)/=0) call mag_turnof

  end select

end subroutine Magnet
      
