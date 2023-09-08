!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_matrsw.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Compute elemental matrix and RHS for the water saturation
!> @details Compute elemental matrix and RHS for the water saturation
!> @} 
!------------------------------------------------------------------------
subroutine por_matrsw()
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
     call por_elmops(1_ip) 
     !
     ! Boundary assembly
     !
     ! call por_bouope()  for the moment we do not need any boundary terms
     !
     ! Add well terms nodally - modifies rhsid for saturation
     !
     call por_wellte()    
  end if

end subroutine por_matrsw
