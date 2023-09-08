subroutine wav_outinf
  !-----------------------------------------------------------------------
  !****f* wavequ/wav_outinf
  ! NAME 
  !    wav_outinf
  ! DESCRIPTION
  !    This routine writes info in the wave equation files
  ! USES
  ! USED BY
  !    wav_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_wavequ
  implicit none

  if(kfl_paral<=0) then
     !
     ! Write information in Result file
     !
     if(kfl_rstar/=2) then

     end if
  end if

end subroutine wav_outinf
      
