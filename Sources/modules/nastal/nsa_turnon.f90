subroutine nsa_turnon
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_turnon
  ! NAME 
  !    nsa_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read some data for the run
  !    - Write some info
  !    - Allocate memory
  ! USES
  !    nsa_openfi
  !    nsa_reaphy
  !    nsa_reabcs
  !    nsa_reanut
  !    nsa_reaous
  !    nsa_outinf
  !    nsa_memall
  ! USED BY
  !    Nastal
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_nastal
  use def_domain
  use mod_nsa_euler, only: EULER, nsa_euler_allocate
  implicit none
  !
  ! Open files
  !
  call nsa_openfi(one)
  !
  ! Initial variables
  !
  call nsa_inivar(zero)
  !
  ! Read the physical problem
  !
  call nsa_reaphy

  !
  ! Read the numerical treatment
  !
  call nsa_reanut

  !
  ! Read the output strategy
  !
  call nsa_reaous
  !
  ! Service: Parall
  !
!!!  call nsa_parall(one)
  !
  ! Initialize some (local) variables
  !
!!!  call nsa_inivar(one)
  !
  ! Read initial and boundary conditions
  !
  call nsa_reabcs
  !
  ! Service: Parall
  !
  call nsa_parall(one)
!  call nsa_parall(two)
  !
  ! Initialize some (local) variables
  !
  call nsa_inivar(one)

  !
  ! Service: Parall
  !
!!!!  call nsa_parall(two)
  !
  ! Initial boundary conditions
  !
  call nsa_inibcs()
  !
  ! Write info
  !
  call nsa_outinf
  !
  ! Allocate work space
  !
  call nsa_memall
  !
  ! Initialize final local variables
  !
  call nsa_inivar(two)
  !
  ! Compute some special stuff for the module
  !
!!!!  call nsa_modspe
  !
  ! Open files needed occasionally
  !
  call nsa_openfi(two)
  !
  ! Close input data file and print some general info on screen
  !
  call nsa_outbcn

  if(euler_nsa) call nsa_euler_allocate( EULER )

end subroutine nsa_turnon
