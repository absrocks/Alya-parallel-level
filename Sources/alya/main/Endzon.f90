!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    Endzon.f90
!> @author  Guillaume Houzeaux
!> @date    30/09/2014
!> @brief   Check block convergence if coupling
!> @details Check block convergence if coupling
!> @} 
!-----------------------------------------------------------------------

subroutine Endzon(itask)
  
  use def_kintyp,    only : ip
  use def_master,    only : iblok
  use def_master,    only : kfl_gocou
  use def_master,    only : mmodu
  use def_master,    only : lmord
  use def_master,    only : itinn
  use def_master,    only : ITASK_ENDZON
  use def_coupli,    only : mcoup
  use def_coupli,    only : kfl_gozon
  use def_coupli,    only : coupling_driver_iteration
  use def_coupli,    only : coupling_driver_number_couplings
  use mod_couplings, only : COU_CHECK_CONVERGENCE
  use mod_moduls,    only : moduls 
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: imodu,iorde
  !
  ! Check convergence if we have coupling
  !
#ifndef COMMDOM
  !
  if( mcoup > 0 ) then
     !
     ! Call modules to exchange data
     !
     call moduls(ITASK_ENDZON)
     !
     ! We are currently in block IBLOK
     !  
     if( coupling_driver_number_couplings(iblok) /= 0 ) then
        call COU_CHECK_CONVERGENCE(iblok,kfl_gozon)
        if( kfl_gozon == 1 ) kfl_gocou = 1
     else
        kfl_gozon = 0 
     end if
     !
     ! Output convergence
     !
     call cou_cvgunk()   
     
  else
     !
     ! Go to next block
     !
     kfl_gozon = 0
  end if
  !
#else 
  call moduls(ITASK_ENDZON)
#if   COMMDOM==1
  kfl_gozon = 0
#endif 
#endif 

end subroutine Endzon
