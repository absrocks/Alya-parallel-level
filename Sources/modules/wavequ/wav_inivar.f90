subroutine wav_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_inivar
  ! NAME 
  !    wav_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    wav_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_wavequ
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( 1_ip )
     !
     ! Postprocess Variable names and types
     !
     postp(1) % wopos(1,1) = 'WAVAM'
     postp(1) % wopos(1,2) = 'WAVVE'
     postp(1) % wopos(1,3) = 'WAVAC'
     postp(1) % wopos(2,1) = 'SCALA'
     postp(1) % wopos(2,2) = 'SCALA'
     postp(1) % wopos(2,3) = 'SCALA'
     cpuit_wav             =  0.0_rp       
     !
     ! Solver
     !     
     call soldef(-1_ip)
     solve(1) % kfl_solve = 1     
     solve(1) % wprob     = 'WAVE_AMPLITUDE'

  case ( 2_ip )
     !
     ! Number of unknown components
     !
     if( kfl_tisch_wav == 2 ) then
        !
        ! Leap-frog scheme scheme
        !
        ncomp_wav = 4

     else if( kfl_tisch_wav == 3 ) then
        !
        ! Newmark scheme
        !
        ncomp_wav = 3

     end if

  end select

end subroutine wav_inivar
