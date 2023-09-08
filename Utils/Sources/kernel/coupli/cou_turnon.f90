!------------------------------------------------------------------------
!> @addtogroup Coupling 
!> @{
!> @file    cou_turnon.f90
!> @date    03/03/2014
!> @author  Guillaume Houzeaux
!> @brief   Turn on coupling 
!> @details Read data and allocate memory
!> @} 
!------------------------------------------------------------------------
subroutine cou_turnon()
  use def_kintyp, only : ip
  implicit none
  !
  ! Initialize variables
  !
  call cou_inivar(1_ip)
  !
  ! Open files
  !  
  call cou_openfi()
  !
  ! Read data
  !  
  call cou_readat()
  !
  ! Broadcast data
  !  
  call cou_parall()
  !
  ! Initialize variables after reading data
  !
  call cou_inivar(2_ip)

end subroutine cou_turnon
