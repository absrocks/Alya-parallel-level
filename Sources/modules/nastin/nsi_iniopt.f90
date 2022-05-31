subroutine nsi_iniopt()
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the adjoint equations of the NS equation.
  !
  !-----------------------------------------------------------------------
  use def_domain
  use def_nastin
  use def_master
  use def_kermod, only : kfl_adj_prob,kfl_cos_opt,sens,sens_mesh,costf,kfl_dvar_type,kfl_cost_type
  use def_solver
  implicit none

  integer(ip)  :: izrhs,idesvar
    
  if (kfl_cos_opt == 1) costf = 0.0_rp
  
  if (kfl_adj_prob == 1) then
!     resdiff_nsi = 0.0_rp 
    dcost_dx_nsi = 0.0_rp
  endif
  
  if (kfl_dvar_type == 5) then
    sens_mesh = 0.0_rp
  endif
  
!   if ( kfl_coupl(ID_NASTIN,ID_TURBUL) == 0 ) then
!     if (kfl_dvar_type == 5) then
!       sens_mesh = 0.0_rp
!     else
!       sens = 0.0_rp
!     endif
!   endif
 
end subroutine nsi_iniopt
 
