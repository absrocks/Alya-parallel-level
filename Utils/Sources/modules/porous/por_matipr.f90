!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_matipr.f90
!> @date    13/08/2013
!> @author  Herbert Owen
!> @brief   Compute elemental matrix and RHS for the pressure inicialization - very similar por_solite but only lapl & grav.
!> @details Compute elemental matrix and RHS for the pressure inicialization - very similar por_solite but only lapl & grav.
!> @} 
!------------------------------------------------------------------------
subroutine por_matipr()
  use def_parame
  use def_master
  use def_porous
  use def_domain
  implicit none
  !
  ! Initializations
  !
  resgs_por(1) = 0.0_rp
  resgs_por(2) = 0.0_rp

  if( INOTMASTER ) then 
     !
     ! Element assembly
     !
     call por_elmop0(1_ip)     
  end if

end subroutine por_matipr
