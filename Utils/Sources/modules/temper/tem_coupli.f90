!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_coupli.f90
!> @author  Guillaume Houzeaux
!> @brief   Coupling of tmper with other modules
!> @details Coupling of tmper with other modules\n
!>          - Coupling with CHEMIC: water vapor model\n
!>            qL = w/Area * Dh \n
!>            w/Area = rho_wv k_wv * gradc_wv.n \n
!> @} 
!-----------------------------------------------------------------------
subroutine tem_coupli(itask)
  use def_master
  use def_domain
  use def_elmtyp
  use def_temper
  use mod_gradie
  use mod_tem_commdom, only: tem_commdom_lm2_code_i, tem_commdom_lm2_code_j
  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( ITASK_CONCOU )
     call tem_commdom_lm2_code_i(-1_ip, -1_ip)
     call tem_commdom_lm2_code_j(-1_ip, -1_ip)

  case ( ITASK_BEGSTE )
     call tem_commdom_lm2_code_i(-1_ip, -1_ip)
     call tem_commdom_lm2_code_j(-1_ip, -1_ip)

  case ( ITASK_INIUNK )
     call tem_commdom_lm2_code_i(-1_ip, -1_ip)
     call tem_commdom_lm2_code_j(-1_ip, -1_ip)

  case ( ITASK_BEGITE )

     if( kfl_coupl(ID_TEMPER,ID_CHEMIC) ==2 ) then
        !
        ! Water vapor concentration gradients
        !
        if( INOTMASTER ) then
           call gradie(conce(1:npoin,1,1),gradc_tem)
        end if

     end if

  case(ITASK_DOITER) 


  end select

end subroutine tem_coupli
