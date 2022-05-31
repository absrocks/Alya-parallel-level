subroutine nsi_begite()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_begite
  ! NAME 
  !    nsi_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the incompressible NS
  !    equations. 
  ! USES
  !    nsi_tittim
  !    sni_updbcs
  !    nsi_inisol
  !    nsi_updunk
  ! USED BY
  !    nsi_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_nsi_hydrostatic,  only : nsi_hydrostatic_pressure
  use mod_messages, only : livinf
  implicit none
  !
  ! Initializations
  !
  kfl_goite_nsi = 1
  itinn(modul)  = 0
  if( momod(modul) % miinn == 0 ) kfl_goite_nsi = 0
  if( itcou == 1 ) call nsi_tistep()
  !
  ! CFD wake
  !
  call nsi_cfdwak(2_ip)
  !
  ! Coupling
  !
  call nsi_coupli(ITASK_BEGITE)
  !
  ! Set up the solver parameters for the NS equations
  !
  call nsi_inisol(one)
  !
  ! Set up the parameters for the optimization
  !
  call nsi_iniopt()
  !
  ! If hydrostatic pressure should be updated
  !
  call nsi_hydrostatic_pressure(ITASK_BEGITE)
  ! 
  ! Force Dirichlet boundary conditions on VELOC(:,:,1)
  !
  call nsi_updbcs(ITASK_BEGITE)
  !
  ! Obtain the initial guess for inner iterations in global
  ! VELOC(:,:,2) <= VELOC(:,:,1)
  !
  call nsi_updunk(ITASK_BEGITE)

end subroutine nsi_begite



