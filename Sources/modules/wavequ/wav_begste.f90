subroutine wav_begste
!-----------------------------------------------------------------------
!****f* Wavequ/wav_begste
! NAME 
!    wav_begste
! DESCRIPTION
!    This routine prepares for a new time step of the temperature
!    equation      
! USES
!    wav_iniunk
!    wav_updtss
!    wav_updbcs
!    wav_updunk
!    wav_radvuf
! USED BY
!    Temper
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  implicit none
  integer(ip), save :: ipass=0

  if(ipass==0) then
     !
     ! First time we pass here: do not use ittim as this can be 
     ! a restart run 
     !
     ipass=1
     call wav_inivar(2_ip)
  end if
 
  if(kfl_stead_wav/=1) then     
     !
     ! Initial guess fo the temperature: u(n,0,*) <-- u(n-1,*,*).
     !
     call wav_updunk(1_ip)
    
  end if
  
end subroutine wav_begste

