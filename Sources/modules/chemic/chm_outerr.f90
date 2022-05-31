subroutine chm_outerr()
!------------------------------------------------------------------------
!****f* Chemic/chm_outerr
! NAME 
!    chm_outerr
! DESCRIPTION
!    This routine checks if there are erros and warnings
! USES
! USED BY
!    chm_turnon
!***
!------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_chemic
  use mod_outfor, only : outfor
  implicit none
  integer(ip)   :: ierro,iwarn
  integer(ip)   :: iclas
  integer(ip)   :: cmean_present
  integer(ip)   :: cvar_present
  integer(ip)   :: chist_present
  integer(ip)   :: zmean_present
  integer(ip)   :: zvar_present
  integer(ip)   :: imean_present

  ierro = 0_ip
  iwarn = 0_ip
  !
  ! TRANSIENT PROBLEM 
  !
  if( kfl_timei /= 0 .and. kfl_timei_chm == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'STEADY CHEMIC IN A TRANSIENT CALCULATION')
  end if
  !
  ! Table framework
  !
  if (kfl_tab_fw_chm > -1) then
     cmean_present = 0_ip
     cvar_present  = 0_ip
     chist_present = 0_ip
     zmean_present = 0_ip
     zvar_present  = 0_ip
     imean_present = 0_ip
     do iclas = 1, table_fw % main_table % ndim
        select case (table_fw % main_table % coords(iclas) % name)
        case ('CMEAN','C    ')
            cmean_present  = iclas
        case ('CVAR ')
            cvar_present   = iclas
        case ('CHIST')
            chist_present  = iclas
        case ('ZMEAN','Z    ')
            zmean_present  = iclas
        case ('ZVAR ')
            zvar_present   = iclas
        case ('IMEAN','I    ')
            imean_present  = iclas
        end select
     enddo

     if(kfl_heat_loss_chm /= 0 .and. imean_present==0 ) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,'I OR IMEAN SHOULD BE IN TABLE FOR NON-ADIABATIC CACULATION')
     end if

     if(kfl_cotur_chm ==0 ) then
        if (cvar_present > 0) then
           if ( table_fw % kfl_scale(cvar_present) /= -1 ) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'PLEASE TURN OFF CVAR SCALING IN TABLE FRAMEWORK')
           endif
        endif
        if (zvar_present > 0) then
           if ( table_fw % kfl_scale(zvar_present) /= -1 ) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu,'PLEASE TURN OFF ZVAR SCALING IN TABLE FRAMEWORK')
           endif
        endif
     end if

     if( kfl_ufpv_chm /= 0 .and. chist_present==0) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,'CHIST SHOULD BE IN TABLE FOR UFPV CACULATION')
     end if

     if( kfl_varZ_chm /= 0 ) then
        if (zvar_present==0) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,'ZVAR SHOULD BE IN TABLE')
        else
            if ( abs(kfl_varZ_chm) == 1 .and.  table_fw % kfl_scale(zvar_present) /= (100+zmean_present) ) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'ZVAR SCALING SHOULD BE VARIANCE IN TABLE FRAMEWORK')
            endif
            if ( abs(kfl_varZ_chm) == 2 .and.  table_fw % kfl_scale(zvar_present) /= -1*(100+zmean_present) ) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'ZVAR SCALING SHOULD BE SQUARE IN TABLE FRAMEWORK')
            endif
        end if
     end if

     if( kfl_varYc_chm /= 0 ) then
        if (cvar_present==0) then
           ierro = ierro + 1
           call outfor(1_ip,momod(modul)%lun_outpu,'CVAR SHOULD BE IN TABLE')
        else
            if (kfl_varYc_chm==1 .and. table_fw % kfl_scale(cvar_present) /= (100+cmean_present) ) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'CVAR SCALING SHOULD BE VARIANCE IN TABLE FRAMEWORK')
            endif
            if (kfl_varYc_chm==2 .and. table_fw % kfl_scale(cvar_present) /= -1*(100+cmean_present) ) then
               ierro = ierro + 1
               call outfor(1_ip,momod(modul)%lun_outpu,'CVAR SCALING SHOULD BE SQUARE IN TABLE FRAMEWORK')
            endif
        end if
     end if

     if( kfl_premix_chm /= 0 .and. zmean_present > 0) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,'ZMEAN SHOULD NOT BE PRESENT FOR PREMIXED CALCULATION')
     end if

     if( kfl_premix_chm == 0 .and. zmean_present == 0) then
        ierro = ierro + 1
        call outfor(1_ip,momod(modul)%lun_outpu,'Z OR ZMEAN SHOULD BE PRESENT FOR NON-PREMIXED CALCULATION')
     end if
  endif

  
  
  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------
  call errors(3_ip,ierro,iwarn,'NULL')

end subroutine chm_outerr
