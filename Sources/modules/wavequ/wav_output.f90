subroutine wav_output()
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_output
  ! NAME 
  !    wav_output
  ! DESCRIPTION
  !    Output of the solution
  ! USED BY
  !    Wavequ
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  use mod_postpr
  use mod_iofile
  implicit none
  integer(ip) :: ivari,ivarp
  !
  ! Initial solution, end of a time step and and of run
  !
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     call wav_outvar(ivari)
  end do

end subroutine wav_output
