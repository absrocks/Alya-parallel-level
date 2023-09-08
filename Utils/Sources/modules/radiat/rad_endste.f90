subroutine rad_endste()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_endste
  ! NAME 
  !    rad_endste
  ! DESCRIPTION
  !    This routine ends a time step of the radiation equation.
  ! USES
  !    rad_cvgunk
  !    rad_updunk
  !    rad_output
  ! USED BY
  !    Radiat
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_radiat
  implicit none
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  !
  ! Update heat radiation source term
  !
  call rad_updunk(7_ip)
  !
  ! Write restart file
  !
  call rad_restar(2_ip)
  !
  ! If not steady, go on
  !
!!F if(kfl_stead_rad==0.and.kfl_timei_rad==1.and.kfl_conve(modul)==1) kfl_gotim = 1  !!F What is this

end subroutine rad_endste
