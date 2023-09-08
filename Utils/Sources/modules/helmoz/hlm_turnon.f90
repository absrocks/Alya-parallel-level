subroutine hlm_turnon()

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_turnon.f90
  ! NAME 
  !    hlm_turnon
  ! DESCRIPTION
  !    This routine starts the run of the 'helmoz' module.
  ! USES
  ! USED BY
  !    Helmoz
  !-----------------------------------------------------------------------

  use def_parame
  use def_master
  use def_helmoz
  implicit none
  !
  ! Initial variables
  !
  call hlm_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call hlm_reaphy()
  !
  ! Initial variables depending on physical problem
  !
  call hlm_inivar(1_ip)
  !
  ! Read the numerical problem
  !
  call hlm_reanut()
  !
  ! Read the postprocess
  !
  call hlm_reaous()
  !
  ! Read the boundary conditions
  !
  call hlm_reabcs()
  !
  ! Parallelization
  !
  call hlm_parall(1_ip)
  !
  ! Impose boundary conditions
  !
  call hlm_inibcs()
  !
  ! Errors and warnings
  !
  call hlm_outerr()
  !
  ! If deflated preconditioning is used, compute the groups that it needs
  !
  if (solve_sol(1) % kfl_preco == 14_ip) call hlm_inivar(2_ip)
  !
  ! Memory
  !
  call hlm_memall()
  !
  ! Compute some arrays
  !
  call hlm_inivar(3_ip)

end subroutine hlm_turnon
