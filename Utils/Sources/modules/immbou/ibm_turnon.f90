subroutine ibm_turnon()
  !-----------------------------------------------------------------------
  !****f* Temper/ibm_turnon
  ! NAME 
  !    ibm_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the temperature equation.
  !    - Write some info
  !    - Check if there are errrors
  !    - Allocate memory
  ! USES
  !    ibm_openfi
  !    ibm_reaphy
  !    ibm_reabcs
  !    ibm_reanut
  !    ibm_reaous
  !    ibm_outinf
  !    ibm_memall
  !    ibm_restar
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  implicit none
  !
  ! Initial variables
  !
  call ibm_inivar(0_ip)
  !
  ! Read the dimension
  !
  call ibm_readim()
  !
  ! Read the geometry
  !
  call ibm_reageo()
  !
  ! IB and Integration
  !
  call ibm_cderda()
  !
  ! Read the physical problem
  !
  call ibm_reaphy()
  !
  ! Read the numerical treatment
  !
  call ibm_reanut()
  !
  ! Read the output and postprocess
  !
  call ibm_reaous()
  !
  ! Parall
  !
  call ibm_parall(1_ip)
  !
  ! Shape functions and derivatives
  !
  call ibm_cshder(1_ip)
  !
  ! Open additional files
  !
  call ibm_openfi(0_ip)

end subroutine ibm_turnon
