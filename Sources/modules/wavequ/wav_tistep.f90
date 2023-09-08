subroutine wav_tistep
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_tistep
  ! NAME 
  !    wav_tittim
  ! DESCRIPTION
  !    This routine sets the time step
  ! USES
  ! USED BY
  !    wav_begite
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_wavequ
  use mod_outfor, only : outfor
  implicit none

  if(kfl_timco/=2) then
     dtinv_wav=dtinv
     if(kfl_stead_wav==1) dtinv_wav = 0.0_rp

     if(kfl_tisch_wav==1) then
        !
        ! Trapezoidal rule: Euler iterations
        !
        if(ittim<=neule_wav) then
           kfl_tiacc_wav=1
        else
           kfl_tiacc_wav=kfl_tiaor_wav
        end if
        if(kfl_tiacc_wav==2) dtinv_wav = 2.0_rp*dtinv_wav
     else
        !
        ! BDF scheme: increase integration order at each time step
        !
        kfl_tiacc_wav=min(kfl_tiaor_wav,ittim)
        call parbdf(kfl_tiacc_wav,pabdf_wav)
     end if
  end if

  routp(1)=dtcri_wav
  routp(2)=0.0_rp
  routp(3)=0.0_rp
  ioutp(1)=1
  ioutp(2)=kfl_stead_wav
  call outfor(8_ip,lun_outpu,' ')

end subroutine wav_tistep
