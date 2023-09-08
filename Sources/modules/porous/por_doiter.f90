!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_doiter.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Controls the internal loop of the porous equations
!> @details Controls the internal loop of the porous equations
!> @} 
!------------------------------------------------------------------------
subroutine por_doiter()
  use def_parame
  use def_master
  use def_solver
  use def_porous
  implicit none

  if( kfl_stead_por == 0 ) then
     call por_begite()
     do while( kfl_goite_por == 1 )
        kprsa_por = 2_ip
        call por_solitw()
        if ( ittim >= ndels_por+1_ip ) call por_endite(three) 
        kprsa_por = 1_ip  
        call por_solite()
        call por_endite(one)
        call por_velsmo()            ! velocity smoothing
     end do
     call por_endite(two)
  end if

end subroutine por_doiter
