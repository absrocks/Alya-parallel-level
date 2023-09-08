
!-----------------------------------------------------------------------
!> @defgroup Exmedi
!> @{
!> @file    Exmedi.f90
!> @author  Mariani Vazques
!> @date    2018-12-30
!> @brief   Main Exmedi surboutine
!> @details This is the main subroutine for the Excitable Media module.
!>          A Bidomain Model is used for cardiac electrical propagation.
!>          Two different kinds of models are implemented for ionic currents:
!>          - No subcell models (approximate models): FitzHugh-Nagumo (FHN)
!>          - Subcell models: TenTusscher (TT), Luo-Rudy (LR), ...
!> @} 
!-----------------------------------------------------------------------

subroutine Exmedi(order)

  use      def_parame
  use      def_domain
  use      def_master
  use      def_solver
  use      def_exmedi

  implicit none
  integer(ip) :: order
  
  select case (order)
     
  case(ITASK_TURNON)
     call exm_turnon
  case(ITASK_INIUNK)
     call exm_iniunk
  case(ITASK_OUTPUT)
     call exm_output
  case(ITASK_TIMSTE) 
     call exm_timste
  case(ITASK_BEGSTE) 
     call exm_begste
  case(ITASK_DOITER)
     call exm_doiter
  case(ITASK_CONCOU)
     call exm_concou
  case(ITASK_CONBLK)
     return
!!     call exm_conblk   --> to be done when needed
  case(ITASK_ENDSTE)
     call exm_endste
  case(ITASK_TURNOF)
     call exm_turnof
     
  end select
  !
  ! Coupling
  ! 
  if( order > 1000 ) call exm_plugin(order-1000_ip) 
  
end subroutine Exmedi
      
