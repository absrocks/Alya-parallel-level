!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_endite.f90
!> @date    29/03/2016
!> @author  Guillaume Houzeaux
!> @brief   End of inner iteration
!> @details End of inner iteration
!> @}
!------------------------------------------------------------------------

subroutine neu_endite(itask)

  use def_parame
  use def_master
  use def_neutro
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask ) 
     
  case ( 1_ip )
     !
     ! Computer residual of current energy and direction
     !
     call neu_cvgunk(1_ip)  !> Compute L2 residual
     call neu_updunk(3_ip)  !> Assign u(n,i,j-1) <-- u(n,i,j), update of unknown after solver
 
  case ( 2_ip )
     !
     ! Compute convergence residual of the external iteration
     !
     call livinf(16_ip,' ',itinn(modul))
     call neu_cvgunk(2_ip)  !> Check convergence of the outer iterations: || T(n,i,*) - T(n,i-1,*)|| / ||T(n,i,*)||
     call neu_updunk(4_ip)  !> Assign u(n,i-1,*) <-- u(n,i,*)

  case ( 3_ip )
     !
     ! Check global convergence over all energies and directions
     !
     call neu_cvgunk(4_ip)  

  end select

end subroutine neu_endite
