!-----------------------------------------------------------------------
!> @defgroup Helmoz
!> @{
!> @file    Helmoz.f90
!> @author  houzeaux
!> @date    2019-01-04
!> @brief   Helmoz
!> @details This routine deals with the incompletely gauged coupled 
!>          vector-scalar potential formulation of Maxwell's equations. 
!>          The task done corresponds to the order given by the master.
!> @} 
!-----------------------------------------------------------------------

subroutine Helmoz(order)

  use def_master

  implicit none

  integer(ip), intent(in) :: order

  select case (order)
    case(ITASK_TURNON)
      call hlm_turnon()
    case(ITASK_TIMSTE) 
      call hlm_timste()
    case(ITASK_INIUNK) 
      call hlm_iniunk()
    case(ITASK_BEGSTE)
      call hlm_begste()
    case(ITASK_DOITER)
      call hlm_doiter()
    case(ITASK_CONCOU)
      call hlm_concou()
    case(ITASK_CONBLK)
      call hlm_conblk()
    case(ITASK_NEWMSH)
      call hlm_newmsh()
    case(ITASK_ENDSTE)
      call hlm_endste()
    case(ITASK_OUTPUT)
      call hlm_output()
    case(ITASK_TURNOF)
      call hlm_turnof()
  end select

end subroutine Helmoz
