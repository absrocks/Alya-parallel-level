subroutine nsi_endste()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_endste
  ! NAME 
  !    nsi_endste
  ! DESCRIPTION
  !    This routine ends a time step of the incompressible NS equations.
  ! USES
  !    nsi_cvgunk
  !    nsi_updunk
  !    nsi_output
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_nastin
  use def_kermod, only : kfl_cos_opt,kfl_adj_prob

  implicit none
  !
  ! Compute convergence residual of the time evolution (that is,
  ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||) and update unknowns
  ! u(n-1,*,*) <-- u(n,*,*) 
  !
  if( kfl_stead_nsi == 0 .and. kfl_timei_nsi == 1 ) then
     call nsi_cvgunk(3_ip) ! Convergence
     call nsi_updunk(5_ip) ! VELOC, PRESS
  end if
  !
  ! Compute averaged variables
  !
  call nsi_averag()
  !
  ! Write restart file
  !
  call nsi_restar(2_ip)
  !
  ! If not steady, go on
  !
  if(kfl_stead_nsi==0.and.kfl_timei_nsi==1.and.kfl_conve(modul)==1) kfl_gotim = 1
  !
  ! Calculate functional and sensitivities
  !
  if (kfl_cos_opt == 1) call nsi_costcal()
  if (kfl_adj_prob == 1) call nsi_senscal()

end subroutine nsi_endste
