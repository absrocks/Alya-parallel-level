subroutine qua_begste
!-----------------------------------------------------------------------
!****f* Quanty/qua_begste
! NAME 
!    qua_begste
! DESCRIPTION
!    This routine prepares for a new time step of the equation      
! USES
!    qua_inivar
!    qua_updbcs
!    qua_updunk
!    qua_dynco
! USED BY
!    Quanty
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  implicit none
  integer(ip), save :: ipass=0

  if(ipass==0) then
     !
     ! First time we pass here: do not use ittim as this can be 
     ! a restart run 
     !
     ipass=1
     ! call qua_inivar(two) ! se usa en restart!!
  end if
 
  if(kfl_stead_qua/=1) then     
     !
     ! Initial guess fo the solution 
     !
     call qua_updunk(one)
     !
     ! Update boundary conditions if it is necesary
     !
     !   call qua_updbcs(one)
     !
     ! Coupling with dynamic solver
     !
     !call qua_dyncou(1_ip)
    
  end if
  
end subroutine qua_begste

