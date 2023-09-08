!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_parall.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Bridge to Parall service
!> @details Bridge to Parall service
!> @} 
!------------------------------------------------------------------------
subroutine por_parall(itask)
  use      def_parame
  use      def_master
  use      def_domain
  use      def_porous
  implicit none
  integer(ip), intent(in) :: itask

  if(kfl_paral>=0) then

     select case(itask)

     case(1)
        !
        ! Exchange data read in por_reaphy, por_reanut and por_reaous
        ! always using MPI, even if this is a partition restart
        !
        call por_sendat(1_ip)

     case(2)
        !
        ! Exchange data read in por_reabcs
        !
!        call por_sendat(2_ip)

     end select

  end if

end subroutine por_parall
