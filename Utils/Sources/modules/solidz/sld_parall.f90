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
  use def_master, only : IPARALL, ISLAVE, NPOIN_INTE_1DIM
  use def_master, only : pari1
  use def_domain, only : npoin
  use def_solidz, only : kacti_sld

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

     case ( 4_ip )
        !
        ! Sum up kacti activation vector of slave neighbors
        !
        if ( ISLAVE ) then
           call vocabu(NPOIN_INTE_1DIM,0_ip,0_ip)
           pari1 => kacti_sld(1:npoin)
           call Parall(400_ip)
        end if

     end select

  end if

end subroutine sld_parall
