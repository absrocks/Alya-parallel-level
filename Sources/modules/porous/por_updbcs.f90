!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_updbcs.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Performs several types of updates for the porous boundary conditions
!> @details Performs several types of updates for the porous boundary conditions
!> @} 
!------------------------------------------------------------------------
subroutine por_updbcs(itask)
  use      def_master
  use      def_domain
  use      def_porous
  implicit none
  integer(ip), intent(in) :: itask
  real(rp), external      :: funcre

  if(kfl_paral/=0) then

     select case(itask)

     case(1)

        !  
        ! Before a time step
        !     

     case(2)
        !
        ! Before a global iteration
        !  
 

     case (3)
        !
        ! Before an inner iteration
        ! 

     end select

  end if

end subroutine por_updbcs
