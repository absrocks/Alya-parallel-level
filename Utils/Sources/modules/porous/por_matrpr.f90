!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_matrpr.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Compute elemental matrix and RHS for the pressure
!> @details Compute elemental matrix and RHS for the pressure
!> @} 
!------------------------------------------------------------------------
subroutine por_matrpr()
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
     call por_elmope(1_ip) 
     !
     ! Boundary assembly
     !
     call por_bouope()  ! boundary terms for gravity
     !
     ! Add well terms nodally - modifies rhsid and amatr for pressure
     !
     call por_wellte() 
    
  end if

end subroutine por_matrpr
