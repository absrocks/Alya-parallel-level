!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_clippi.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   Clipping
!> @details Clipping
!> @}
!-----------------------------------------------------------------------

subroutine chm_clippi()

  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_memory,      only     : memory_alloca, memory_deallo
  use mod_interp_tab,  only     : fw_scale_cont_var 
  use def_chemic,      only     : table_fw
  implicit none
  integer(ip)               :: ipoin,kpoin,iclas
  real(rp) , pointer        :: control(:)   ! input of table lookup function 
  real(rp) , pointer        :: scale_control(:)
  real(rp) , pointer        :: lim_control(:,:)

  select case(kfl_model_chm)

  case (1) ! At Begste

     !
     ! Prevent unde/over shoots
     !
     if ( kfl_negat_chm == 1 .or. kfl_posit_chm == 1) then

        !
        ! Unreacting spray
        !
        if (kfl_premix_chm == 1 .and. kfl_spray_chm /= 0_ip ) then

           kpoin = 0
           do ipoin = 1,npoin
              !
              ! Liquid volume fraction: phi_L
              !
              kpoin = (ipoin - 1) * nclas_chm + 3
              unkno(kpoin) = max(0.0_rp,min(1.0_rp,unkno(kpoin)))
           end do

        else

           nullify(control)
           nullify(scale_control)
           nullify(lim_control)
           call memory_alloca(mem_modul(1:2,modul),'control',      'chm_outvar',control,      5_ip)
           call memory_alloca(mem_modul(1:2,modul),'scale_control','chm_outvar',scale_control,5_ip)
           call memory_alloca(mem_modul(1:2,modul),'lim_control',  'chm_outvar',lim_control,  5_ip,2_ip)

           kpoin = 0
           do ipoin = 1,npoin
              kpoin = (ipoin - 1) * nclas_chm 

              do iclas = 1, table_fw % main_table % ndim
                 select case (table_fw % main_table % coords(iclas) % name)
                 case ('CMEAN','C    ')
                     control(iclas) = unkno(kpoin+1_ip)
                 case ('CVAR ')
                     control(iclas) = unkno(kpoin+2_ip)
                 case ('CHIST')
                     control(iclas) = xZr_chm(ipoin) + xZs_chm(ipoin)
                 case ('ZMEAN','Z    ')
                     control(iclas) = unkno(kpoin+3_ip)
                 case ('ZVAR ')
                     control(iclas) = unkno(kpoin+4_ip)
                 case ('IMEAN','I    ')
                     control(iclas) = therm(ipoin,1)
                 end select
              enddo

              call fw_scale_cont_var( control, scale_control, lim_control, table_fw)
 
              do iclas = 1, table_fw % main_table % ndim
                 select case (table_fw % main_table % coords(iclas) % name)
                 case ('CMEAN','C    ')
                     unkno(kpoin+1) = max(lim_control(iclas,1),min(lim_control(iclas,2),unkno(kpoin+1)))
                 case ('CVAR ')
                     unkno(kpoin+2) = max(lim_control(iclas,1),min(lim_control(iclas,2),unkno(kpoin+2)))
                 case ('ZMEAN','Z    ')
                     unkno(kpoin+3) = max(lim_control(iclas,1),min(lim_control(iclas,2),unkno(kpoin+3)))
                 case ('ZVAR ')
                     unkno(kpoin+4) = max(lim_control(iclas,1),min(lim_control(iclas,2),unkno(kpoin+4)))
                 end select
              enddo

           end do

           call memory_deallo(mem_modul(1:2,modul),'control','chm_outvar',control)
           call memory_deallo(mem_modul(1:2,modul),'scale_control','chm_outvar',scale_control)
           call memory_deallo(mem_modul(1:2,modul),'lim_control','chm_outvar',lim_control)

           if (kfl_spray_chm /= 0_ip ) then
              kpoin = 0
              do ipoin = 1,npoin
                 !
                 ! Liquid volume fraction: phi_L
                 !
                 kpoin = (ipoin - 1) * nclas_chm + 5
                 unkno(kpoin) = max(0.0_rp,min(1.0_rp,unkno(kpoin)))
              end do
           end if

        end if

     end if

  end select

end subroutine chm_clippi
