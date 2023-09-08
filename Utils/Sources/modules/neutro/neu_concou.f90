!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_concou.f90
!> @date    01/04/2016
!> @author  Guillaume Houzeaux
!> @brief   Finalize coupling iteration
!> @details Check convergence 
!> @} 
!-----------------------------------------------------------------------

subroutine neu_concou()

  use def_parame
  use def_master
  use def_neutro
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if( kfl_conve(modul) == 1 .and. resid_neu > cotol_neu ) kfl_gocou = 1
  glres(modul) = resid_neu
  !
  ! Output residuals
  !
  coutp(1) = 'RADIATION'
  routp(1) = resid_neu
  call outfor(9_ip,lun_outpu,' ')

end subroutine neu_concou
