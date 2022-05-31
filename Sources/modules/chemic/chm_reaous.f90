subroutine chm_reaous()
  !------------------------------------------------------------------------
  !****f* partis/chm_reaous
  ! NAME 
  !    chm_reaphy
  ! DESCRIPTION
  !    This routine reads the output strategy
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_chemic
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: dummi

  if(kfl_paral<=0) then
     !
     ! Initializations
     !
     ipara_chm     = 0       ! Integer parameters
     rpara_chm     = 0.0_rp  ! Real parameters
     droplet_postprocess_frequency_chm = 1_ip ! Droplet postprocess frequency
     !
     ! Reach the section
     !
     call ecoute('chm_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('chm_reaous')
     end do
     !
     ! Begin to read data.
     !
     do while(words(1)/='ENDOU')
        call ecoute('chm_reaous')
        
        call posdef(2_ip,dummi) 

        if( words(1) == 'SPECY' ) then
           ipara_chm(1)=getint('SPECY',1_ip,'#Specy number')
        else if ( words(1) == 'FREQD' ) then
           droplet_postprocess_frequency_chm = getint('FREQD',1_ip,'#Droplet postprocess frequency')
           if ( droplet_postprocess_frequency_chm == 0 ) then
              call runend('CHEMIC REAOUS: Droplet postrprocess frequency should not be zero, look at your chm.dat file')
           end if
        end if

     end do

  end if

end subroutine chm_reaous
    
