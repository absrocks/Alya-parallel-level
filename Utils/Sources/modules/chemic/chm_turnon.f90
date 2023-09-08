subroutine chm_turnon()
  !-----------------------------------------------------------------------
  !****f* partis/chm_turnon
  ! NAME 
  !    chm_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the ADS equation.
  !    - Write some info
  !    - Check if there are errrors
  !    - Allocate memory
  ! USES
  !    chm_openfi
  !    chm_reaphy
  !    chm_reabcs
  !    chm_reanut
  !    chm_reaous
  !    chm_outinf
  !    chm_memall
  !    chm_restar
  ! USED BY
  !    Partis
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_chemic
  implicit none
  !
  ! Initial variables
  !
  call chm_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call chm_reaphy()
  !
  ! Postprocess variables
  !
  call chm_inivar(1_ip)
  !
  ! Read the numerical treatment
  !
  call chm_reanut()
  !
  ! Read the output strategy
  !
  call chm_reaous()
  !
  ! Parall service
  !
  call chm_parall(1_ip)
  !
  ! Initial variables
  !
  call chm_inivar(3_ip)
  !
  ! Reaction graph
  !
  call chm_grarea()
  !
  ! Check reactions
  !
  call chm_chkrea()
  !
  ! Allocate memory
  !
  call chm_memall()
  !
  ! Read the boundary conditions
  !
  call chm_reabcs()
  !
  ! Parall service
  !
  call chm_parall(2_ip)
  !
  ! Impose boundary conditions
  !
  call chm_inibcs()
  !
  ! Initial variables
  !
  call chm_inivar(2_ip)
  !
  ! Split assembly arrays
  !
  call chm_splarr()
  !
  ! Write info
  !
  call chm_outinf()
  !
  ! Warnings and errors
  !
  call chm_outerr()
  !
  ! Open additional files
  !
  call chm_openfi(2_ip)

end subroutine chm_turnon
