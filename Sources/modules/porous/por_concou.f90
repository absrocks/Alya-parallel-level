!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_concou.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Checks the porous convergence of the run  --- Sill missing to finish
!> @details Checks the porous convergence of the run
!> @} 
!------------------------------------------------------------------------
subroutine por_concou()
  use def_parame
  use def_master
  use def_porous
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if ( (resid_por(1) > cotol_por) .or. (resid_por(2) > cotol_por) ) kfl_gocou = 1
  end if
  glres(modul) = max(resid_por(1),resid_por(2))
  !
  ! Output residuals
  !
  coutp(1) = 'POROUS'
  routp(1) = sqrt(dot_product(resid_por,resid_por))
  call outfor(9_ip,lun_outpu,' ')

end subroutine por_concou
