subroutine tur_iniopt()
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the adjoint equations of the turbulene equation.
  !
  !-----------------------------------------------------------------------
  use def_domain
  use def_turbul
  use def_master
  use def_kermod, only : kfl_dvar_type,sens,sens_mesh,kfl_cos_opt,costf
  use def_solver
  implicit none

  integer(ip)  :: izrhs,idesvar
  
  if (kfl_cos_opt == 1) costf = 0.0_rp
  if (kfl_dvar_type == 6) then
    sens_mesh = 0.0_rp
  endif
 
end subroutine tur_iniopt
