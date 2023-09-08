!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_begste.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Prepares for a new time step
!> @details Prepares for a new time step
!> @} 
!------------------------------------------------------------------------
subroutine por_begste()
  use def_parame
  use def_master
  use def_domain
  use def_porous
  implicit none
  integer(ip), save :: ipass=0

  if(ipass==0) then
     !
     ! First time we pass here: do not use ittim as this can be 
     ! a restart run 
     !
     ipass=1
!     call por_inivar(two)  This is for subscales - not ready
  end if

  if(kfl_stead_por/=1) then     
     !
     ! Initial guess for the press and Water sat: P(n,0,*) <-- P(n-1,*,*).
     !
     call por_updunk(one)
     !
     ! Update boundary conditions
     !
     call por_updbcs(one)

  end if

end subroutine por_begste

