!------------------------------------------------------------------------
!> @addtogroup Porous
!> @{
!> @file    Porous.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   2 Phase compressible porous flow
!> @details 2 Phase compressible porous flow main subroutine
!> @}
!------------------------------------------------------------------------
subroutine Porous(order)
  use def_master
  use def_kintyp
!  use def_solver
  implicit none
  integer(ip), intent(in) :: order

  select case (order)

  case(ITASK_TURNON)
     call por_turnon()
  case(ITASK_TIMSTE) 
     call por_timste()
  case(ITASK_INIUNK) 
     call por_iniunk()
  case(ITASK_BEGSTE) 
     call por_begste()
  case(ITASK_DOITER)
     call por_doiter()
  case(ITASK_CONCOU)
     call por_concou()
  case(ITASK_CONBLK)
     call por_conblk()
  case(ITASK_NEWMSH)
!     call por_newmsh()    ! not ready
  case(ITASK_ENDSTE)
     call por_endste()
  case(ITASK_OUTPUT)
     call por_output()
  case(ITASK_TURNOF)
     call por_turnof()

  end select

end subroutine Porous
