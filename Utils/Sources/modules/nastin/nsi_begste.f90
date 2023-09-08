subroutine nsi_begste()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_begste
  ! NAME 
  !    nsi_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the incompressible NS
  !    equations.
  ! USES
  !    nsi_iniunk
  !    nsi_updtss
  !    nsi_updbcs
  !    nsi_updunk
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin

  implicit none

  if( kfl_stead_nsi /= 1 )  then     
     !
     ! Initial guess fo the velocity: u(n,0,*) <-- u(n-1,*,*).
     !
     call nsi_updunk(1_ip)
     !
     ! Update frame of reference
     !
     call nsi_updfor()
     !
     ! Update boundary conditions
     !
     call nsi_updbcs(1_ip)
     !
     ! Initialize some variables before starting the time step
     !
     call nsi_inivar(3_ip)
     !
     ! Coupling with dynamic solver
     !
     call nsi_dyncou(1_ip)

  end if 
  !
  ! Surface tension
  !
  if( kfl_surte_nsi /= 0 ) call nsi_norcur

  call nsi_coupli(ITASK_BEGSTE)
  !
  ! Initialize coupling
  !
  call nsi_plugin_init()

end subroutine nsi_begste
