subroutine got_endste()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_endste
  ! NAME 
  !    got_endste
  ! DESCRIPTION
  !    This routine ends a time step of the incomcdropible NS equations.
  ! USES
  !    got_cvgunk
  !    got_updunk
  !    got_output
  ! USED BY
  !    Gotita
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_gotita
  implicit none
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if(kfl_stead_got==0.and.kfl_timei_got==1) then
     call got_cvgunk(three)
     call got_updunk(five)
  end if
  !
  ! Write restart file
  !
  call got_restar(two)
  !
  ! If not steady, go on
  !
  if(kfl_stead_got==0) then
     if(kfl_timei_got==1) kfl_gotim = 1
  end if

end subroutine got_endste
