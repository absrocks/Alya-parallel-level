subroutine chm_scale_Yc()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_scale_Yc
  ! NAME 
  !    chm_scale_Yc
  ! DESCRIPTION
  !    Initialization of the reaction progress variable with equilibrium values
  !
  ! USES
  ! USED BY
  !    chm_iniunk
  !   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only          : ip,rp
  use def_master, only          : conce, therm
  use def_domain, only          : npoin
  use def_chemic, only          : kfl_fields_scale_chm,ncomp_chm, &
                                  kfl_heat_loss_chm,xZr_chm,xZs_chm
  use mod_interp_tab, only      : fw_scale_cont_var 
  use def_chemic,     only      : table_fw
  implicit none
  integer(ip)               :: ipoin,iclas,ind_cmean 
  real(rp)                  :: z,c,z_var,chi_st
  real(rp)                  :: y_c_eq,y_c_0

  real(rp)                  :: control(table_fw % main_table % ndim)          
  real(rp)                  :: scale_control(table_fw % main_table % ndim)     
  real(rp)                  :: lim_control(table_fw % main_table % ndim,2_ip)   
  

  !
  ! Initialization given by geometrical fields with c = 1 or c = 0
  !
  do ipoin = 1,npoin

        control = 0.0_rp
        do iclas = 1, table_fw % main_table % ndim
           select case (table_fw % main_table % coords(iclas) % name)
           case ('CMEAN','C    ')
               control(iclas) = conce(ipoin,1,1)
               ind_cmean = iclas
           case ('CVAR ')
               control(iclas) = conce(ipoin,2,1)
           case ('CHIST')
               control(iclas) = xZr_chm(ipoin) + xZs_chm(ipoin)
           case ('ZMEAN','Z    ')
               control(iclas) = conce(ipoin,3,1)
           case ('ZVAR ')
               control(iclas) = conce(ipoin,4,1)
           case ('IMEAN','I    ')
               if (kfl_heat_loss_chm /= 0) control(iclas) = therm(ipoin,1)
           end select
        enddo

        call fw_scale_cont_var( control, scale_control, lim_control, table_fw)

        select case (kfl_fields_scale_chm)
            !
            ! Impose equilibrium value: y_c_eq 
            !
            case (1_ip)
                if ( conce(ipoin,1,1) > 0.0_rp ) then
                    conce(ipoin,1,1:ncomp_chm) = lim_control(ind_cmean,2)
                else
                    conce(ipoin,1,1:ncomp_chm) = 0.0_rp
                end if
            
            !
            ! convert c to Yc
            !
            case (2_ip)
                conce(ipoin,1,1:ncomp_chm) =lim_control(ind_cmean,1) + conce(ipoin,1,1) * (lim_control(ind_cmean,2) - lim_control(ind_cmean,1))
        end select
  end do
end subroutine chm_scale_Yc
