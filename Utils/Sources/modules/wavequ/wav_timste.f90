subroutine wav_timste
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_timste
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
  !
  ! For the first time step, set up the initial condition
  !
  if(ittim==0) then
     kfl_stead_wav = 0
  end if
  !
  ! Time step size 
  !
  if(kfl_stead_wav/=1) call wav_updtss()     

end subroutine wav_timste

