subroutine insitu(order)
!-----------------------------------------------------------------------
!****f* INSITU VIZ
! NAME 
!    insitu
! DESCRIPTION
!    This is the main subroutine for insitu visulization based on NVIDIA
!    INDEX
!
!***
!
!    The task done corresponds to the order given by the master.
!-----------------------------------------------------------------------

  use def_insitu
  use def_master
  use def_domain
  implicit none
  integer(ip) :: order

  select case (order)

  case(ITASK_TURNON)
     call ins_turnon
     
  case(ITASK_BEGSTE)
     call ins_begste
     
  case(ITASK_ENDSTE)
     call ins_endste
     
  case(ITASK_TURNOF)
     call ins_turnof
     
  end select
  !
  ! Coupling
  ! 
  if( order > 1000 ) call ins_plugin(order-1000_ip) 

end subroutine insitu
