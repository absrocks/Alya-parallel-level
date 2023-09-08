subroutine wav_endste
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_endste
  ! NAME 
  !    wav_endste
  ! DESCRIPTION
  !    This routine ends a time step of the temperature equation.
  ! USES
  !    wav_cvgunk
  !    wav_updunk
  !    wav_output
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_wavequ
  implicit none

  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  if(kfl_stead_wav==0) then
     call wav_cvgunk(three)
     call wav_updunk( five)
  end if
  !
  ! Write restart file
  !
  call wav_restar(two)
  !
  ! If not steady, go on
  !
  if(kfl_stead_wav==0) kfl_gotim = 1

end subroutine wav_endste
