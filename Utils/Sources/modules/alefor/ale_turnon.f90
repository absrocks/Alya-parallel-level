subroutine ale_turnon
  !-----------------------------------------------------------------------
  !****f* Temper/ale_turnon
  ! NAME 
  !    ale_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the ALE formulation equation.
  !    - Write some info
  !    - Allocate memory
  ! USES
  !    ale_openfi
  !    ale_reaphy
  !    ale_reabcs
  !    ale_reanut
  !    ale_reaous
  !    ale_outinf
  !    ale_memall
  ! USED BY
  !    Alefor
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_alefor

  implicit none
  !
  ! Initial variables
  !
  call ale_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call ale_reaphy()
  !
  ! Read the numerical treatment
  !
  call ale_reanut()
  !
  ! Read the output strategy
  !
  call ale_reaous()
  !
  ! Read the boundary conditions
  !
  call ale_reabcs()
  !
  ! Parall service
  !
  call ale_parall(1_ip)
  !
  ! Initial variables
  !
  call ale_inibcs()
  !
  ! Initial variables
  !
  call ale_inivar(1_ip)
  !
  ! Allocate memory
  !
  call ale_memall()

end subroutine ale_turnon
