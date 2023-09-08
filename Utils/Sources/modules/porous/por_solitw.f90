!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_solitw.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Solves an iteration of the porous equations - water saturation.
!> @details Solves an iteration of the porous equations - water saturation.
!> @} 
!------------------------------------------------------------------------
subroutine por_solitw()
  use def_parame
  use def_master
  use def_domain
  use def_porous
  use mod_gradie
  implicit none
  !
  ! Update inner iteration counter and write headings in the solver file.
  !
  itinn(modul) = itinn(modul) + 1
  ittot_por    = ittot_por + 1
  if ( ittim >= ndels_por+1_ip ) then
     !
     ! set solver to water saturation solver
     !
     call por_inisol()
     !
     ! Water Saturation - Construct the right-hand-side
     !
     call por_matrsw()   ! similar to por_matrpr but for the moment only to assemble a rhs for explicit solution
     !
     ! initialize unkno to Sw for solver - actually for the explicit case I do not need it
     !
     call por_updunk(8_ip)
     !
     ! Solve the algebraic system - actually it is explicit but it is better to treat it as if it were a solver 
     !                                                       - see  nsa_inivar L172   solve(1) % kfl_algso = 9
     !
     call solver(rhsid,unkno,amatr,pmatr)  ! here actually we could send dummi for amatr since it is not used
  end if

end subroutine por_solitw
