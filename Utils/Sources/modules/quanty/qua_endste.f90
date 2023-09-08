subroutine qua_endste()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_endste
  ! NAME 
  !    qua_endste
  ! DESCRIPTION
  !    This routine ends a time step of the Shrodinger equation.
  ! USES
  !    qua_cvgunk
  !    qua_updunk
  !    qua_output
  !    qua_restar
  !    qua_dyncou
  ! USED BY
  !    Quanty
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_quanty
  implicit none
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if( kfl_stead_qua == 0 .and. kfl_timei_qua == 1 ) then
     call qua_cvgunk(three)
     call qua_updunk( five)
  end if
  !
  ! Write restart file
  !
  ! call qua_restar(two)
  !
  ! Coupling with dynamic solver
  !
  !call qua_dyncou(2_ip)
  !
  ! If not steady, go on
  !
  if(kfl_stead_qua==0) then
     if(kfl_timei_qua==1) kfl_gotim = 1
  end if

end subroutine qua_endste
