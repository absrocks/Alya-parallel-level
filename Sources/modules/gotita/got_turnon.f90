subroutine got_turnon
  !-----------------------------------------------------------------------
  !****f* Gotita/got_turnon
  ! NAME 
  !    got_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the incomcdropible NS equations.
  !    - Write some info
  !    - Check of there are errors
  !    - Allocate memory
  ! USES
  !    got_openfi
  !    got_reaphy
  !    got_reabcs
  !    got_reanut
  !    got_reaous
  !    got_outinf
  !    got_outerr
  !    got_memall
  ! USED BY
  !    Gotita
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip
  use def_master
  implicit none
  !
  ! Initial variables
  !
  call got_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call got_reaphy()
  !
  ! Read the numerical treatment
  !
  call got_reanut()
  !
  ! Read the output strategy
  !
  call got_reaous()
  !
  ! Service: Parall
  !
  call got_parall(1_ip)
  !
  ! Initial variables
  !
  call got_inivar(1_ip)
  !
  ! Read the boundary conditions
  !
  call got_reabcs()
  !
  ! Service: Parall
  !
  call got_parall(2_ip)
  !
  ! Initial variables
  !
  call got_inivar(2_ip)
  !
  ! Write info
  !
  call got_outinf()
  !
  ! Warnings and errors
  !
  call got_outerr()
  !
  ! Allocate memory
  !
  call got_memall()
  !
  ! Read restart file
  !
  call got_restar(1_ip)
  !
  ! Open additional files
  !
  call got_openfi(2_ip)

end subroutine got_turnon
