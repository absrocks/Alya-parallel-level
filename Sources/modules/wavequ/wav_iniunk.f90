subroutine wav_iniunk
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_iniunk
  ! NAME 
  !    wav_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the wave amplitude
  ! USED BY
  !    wav_begste
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_wavequ
  use      mod_memchk
  implicit none
  integer(ip) :: ipoin,istat

  if( kfl_rstar==0 ) then

     if(kfl_paral/=0) then
        !
        ! Load initial conditions for the wave amplitude
        !
        if(kfl_onnod_wav==1) then
           do ipoin=1,npoin
              wavam(ipoin,3)         = bvess_wav(ipoin)
              wavam(ipoin,ncomp_wav) = bvess_wav(ipoin)
              wavam(ipoin,1)         = bvess_wav(ipoin)
           end do
        end if
     end if

  else
     !
     ! restart
     !
     call wav_restar(1_ip)
  end if

end subroutine wav_iniunk
