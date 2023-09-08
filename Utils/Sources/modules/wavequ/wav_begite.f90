subroutine wav_begite
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_begite
  ! NAME 
  !    wav_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the temperature
  !    equation
  ! USES
  !    wav_tittim
  !    wav_updbcs
  !    wav_inisol
  !    wav_updunk
  ! USED BY
  !    wav_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_wavequ
  use mod_messages, only : livinf
  implicit none
  !
  ! Initializations
  !
  kfl_goite_wav = 1 
  itinn(modul)     = 0
  if(itcou==1) call wav_tistep
  call livinf(15_ip,' ',modul)
  !
  ! Set up the solver parameters for the temperature equation
  !
  call wav_inisol  
  !
  ! Obtain the initial guess for inner iterations
  !
  call wav_updunk(2_ip)

end subroutine wav_begite
    
