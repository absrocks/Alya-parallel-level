subroutine Nastal(order)
!-----------------------------------------------------------------------
!****f* Nastal/nastal
! NAME 
!    Nastal
! DESCRIPTION
!    This routine deals with the full set of NS equations. 
! USES
!    nsa_turnon
!    nsa_begste
!    nsa_doiter
!    nsa_gencon
!    nsa_endste
!    nsa_turnof
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
     if(kfl_modul(modul)/=0) call nsa_turnon
  case(ITASK_INIUNK)
     if(kfl_modul(modul)/=0) call nsa_iniunk
  case(ITASK_OUTPUT)
     if(kfl_modul(modul)/=0) call nsa_output
  case(ITASK_TIMSTE) 
     if(kfl_modul(modul)/=0) call nsa_timste
  case(ITASK_BEGSTE)
     if(kfl_modul(modul)/=0) call nsa_begste
  case(ITASK_DOITER)
     if(kfl_modul(modul)/=0) call nsa_doiter
  case(ITASK_CONCOU)
     if(kfl_modul(modul)/=0) call nsa_concou
  case(ITASK_CONBLK)
     if(kfl_modul(modul)/=0) call nsa_conblk
  case(ITASK_NEWMSH)
!     if(kfl_modul(modul)/=0) call nsa_newmsh
  case(ITASK_ENDSTE)
     if(kfl_modul(modul)/=0) call nsa_endste
  case(ITASK_TURNOF)
     if(kfl_modul(modul)/=0) call nsa_turnof

  end select
  !
  ! Coupling
  ! 
  if( order > 1000 ) call nsa_plugin(order-1000_ip) 

end subroutine Nastal
      
