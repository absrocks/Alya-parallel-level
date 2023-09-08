subroutine wav_concou
!-----------------------------------------------------------------------
!****f* Wavequ/wav_concou
! NAME 
!    wav_concou
! DESCRIPTION
!    This routine checks the wave equation convergence of the run.
! USED BY
!    Wavequ
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_wavequ
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if(resid_wav>cotol_wav) kfl_gocou = 1
  end if
  glres(modul) = resid_wav
  !
  ! Output residuals
  !
  coutp(1)='WAVE AMPLITUDE'
  routp(1)=resid_wav
  call outfor(9_ip,lun_outpu,' ')

end subroutine wav_concou
