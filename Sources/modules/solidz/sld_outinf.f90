!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_builtin_materials.f90
!> @author  Solidz team
!> @date    May 2020
!> @brief   
!> @details 
!> @}
!-----------------------------------------------------------------------
subroutine sld_outinf(itask)
  use def_kintyp_basic,      only : ip, rp, lg
  use def_master,            only : INOTSLAVE, kfl_rstar
  use def_solidz,            only : kfl_volca_sld
  use mod_sld_cardiac_cycle, only : sld_cardiac_cycle_write_res
  implicit none
  integer(ip), intent(in) :: itask

  if( INOTSLAVE ) then

    select case(itask)

      case( 1_ip )
        !
        ! Writhing heading files
        !
        if( kfl_rstar /= 2 )then
          if( kfl_volca_sld > 0_ip ) call sld_cardiac_cycle_write_res(1_ip)
        end if

      case( 2_ip )
        !
        ! Write results at the end of the step
        !
        if( kfl_volca_sld > 0_ip ) call sld_cardiac_cycle_write_res(2_ip)

    end select

  end if

end subroutine sld_outinf
      
