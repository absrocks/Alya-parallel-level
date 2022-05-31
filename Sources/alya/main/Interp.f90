!-----------------------------------------------------------------------
!> @addtogroup Interp
!> @{
!> @file    Interp.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Interpribution
!> @details interpribution of data
!> @} 
!-----------------------------------------------------------------------

subroutine Interp()
  
  use def_kintyp,   only : ip
  use def_master,   only : ITASK_INTERP,kfl_paral
  use def_master,   only : iblok
  use def_master,   only : nblok
  use mod_moduls,   only : moduls
  use mod_messages, only : messages_live
  implicit none
  
  call messages_live('INTERPOLATION','START SECTION')
  do iblok = 1,nblok
     call moduls(ITASK_INTERP)
  end do
  call Kermod(ITASK_INTERP)
  call messages_live('INTERPOLATION','END SECTION')
    
end subroutine Interp
