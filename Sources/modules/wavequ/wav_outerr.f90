subroutine wav_outerr
  !------------------------------------------------------------------------
  !****f* Wavequ/wav_outerr
  ! NAME 
  !    wav_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    wav_turnon
  !***
  !------------------------------------------------------------------------
  use      def_master
  use      def_wavequ
  use mod_outfor, only : outfor
  implicit none
  integer(ip)   :: ierro=0,iwarn=0
  character(20) :: messa
  !
  ! Compatibility of time schemes
  !
  if(kfl_timet_wav==1) then
     if(kfl_tisch_wav/=2) then
        ierro=ierro+1
        call outfor(1_ip,momod(modul)%lun_outpu,&
             'ONLY LEAP-FROG SCHEME AVAILABLE IN EXPLICIT METHOD')
     end if
  end if
  !
  ! Warning
  !
  if(iwarn/=0) call outfor(3_ip,momod(modul)%lun_outpu,' ')
  !
  ! Stop
  !
  messa=intost(ierro)
  if(ierro/=0) call outfor(4_ip,momod(modul)%lun_outpu,trim(messa))

end subroutine wav_outerr
