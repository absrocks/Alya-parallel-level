!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_parall.f90
!> @author  Solidz
!> @date
!> @brief   This routine is a bridge to Parall service
!> @details
!> @}
!-----------------------------------------------------------------------

subroutine sld_parall(itask)

  use def_kintyp, only : ip
  use def_master, only : IPARALL

  implicit none

  integer(ip), intent(in) :: itask !< What to do

   if ( IPARALL ) then

      select case ( itask )

      case ( 1_ip )
         !
         ! Exchange data read in sld_reaphy, sld_reanut and sld_reaous
         ! always using MPI, even if this is a partition restart
         !
         call sld_sendat(1_ip)

      end select

   end if

end subroutine sld_parall
