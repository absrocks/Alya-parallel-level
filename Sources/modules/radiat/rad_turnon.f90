subroutine rad_turnon()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_turnon
  ! NAME 
  !    rad_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the P1 equation.
  !    - Write some info
  !    - Check if there are errrors
  !    - Allocate memory
  ! USES
  !    rad_openfi
  !    rad_reaphy
  !    rad_reabcs
  !    rad_reanut
  !    rad_reaous
  !    rad_outinf
  !    rad_memall
  !    rad_restar
  ! USED BY
  !    Radiat
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_radiat
  implicit none
  !
  ! Initial variables
  !
  call rad_inivar(zero) 
  !
  ! Read the physical problem
  !
  call rad_reaphy() 
  !
  ! Read the numerical treatment
  !
  call rad_reanut() !!!F  Still things I don´t know what they are
  !
  ! Read the output strategy
  !
  call rad_reaous()
  !
  ! Read the boundary conditions
  !
  call rad_reabcs() 
  !
  ! Parall service
  !
  call rad_parall(1_ip)
  !
  ! Initial variables
  !
  call rad_inibcs() 
  !
  ! Initial variables
  !
  call rad_inivar(1_ip)
  !
  ! Write info
  !
  call rad_outinf() !!!F
  !
  ! Warnings and errors
  !
  call rad_outerr()  !!!F Add my own error checking
  !
  ! Allocate memory
  !
  call rad_memall() 
  !
  ! Open additional files
  !
  call rad_openfi(two) !!!F Still some files I don´t understand

end subroutine rad_turnon
