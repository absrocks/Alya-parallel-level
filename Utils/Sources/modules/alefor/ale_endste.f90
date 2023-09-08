subroutine ale_endste()
!------------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_endste.f90
!> @author  Guillaume Houzeaux
!> @date    February, 2017
!> @brief   Finish the time step
!>
!> @details Finish the time step
!>
!> @}
!------------------------------------------------------------------------

  use def_parame
  use def_master
  use def_alefor
  implicit none
  !
  ! Update RB unknowns
  !
  call ale_updunk(25_ip)
  !
  ! Write restart file
  !
  call ale_restar(2_ip)
  !
  ! Update unknowns d(n-1,*,*) >-- d(n,*,*)
  !
  call ale_updunk(5_ip)
  !
  ! When there is a timefunction and ALEFOR is running alone, check time advance
  !
  if(kfl_timef_ale==1) kfl_gotim = 1

end subroutine ale_endste
