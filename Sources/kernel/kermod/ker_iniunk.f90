subroutine ker_iniunk()
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_iniunk
  ! NAME 
  !    ker_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the velocity.
  !    If this is a restart, initial condition are loaded somewhere else.
  !    Values are stored in position
  !    VELOC(:,:,NPREV_NSI) 
  !    PRESS(:,:,NPREV_NSI) 
  ! USED BY
  !    ker_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame 
  use def_master
  use def_domain
  use def_kermod
  use mod_ker_proper
  implicit none
  !
  ! Compute wall distance and wall normal
  !
  call ker_walgen(1_ip)
  call ker_walnor(1_ip)
  !
  ! Calculate roughness
  !
  call ker_roughn()
  !
  ! Calculate canopy height & height over terrain
  !
  call ker_canopy()
  !
  ! Velocity and temperature field
  !
  call ker_velfun(ITASK_INIUNK)
  call ker_temfun(ITASK_INIUNK)
  call ker_confun(ITASK_INIUNK)
  call ker_disfun(ITASK_INIUNK)
  ! 
  ! Properties
  !
  call ker_proper_allocate_properties()
  call ker_proper_check_properties()     
  call ker_updpro(ITASK_INIUNK)
  
  !
  ! Test properties
  !
  !call ker_tespro()
!  call ker_waexlo()
  
end subroutine ker_iniunk
