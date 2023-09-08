subroutine dod_turnon
  !-----------------------------------------------------------------------
  !****f* dodeme/dod_turnon
  ! NAME 
  !    dod_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for Dodeme strategy.
  !    - Write some info
  !    - Check of there are errors
  !    - Allocate memory
  ! USES
  !    dod_openfi
  !    dod_readat
  ! USED BY
  !    Dodeme
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  implicit none 
  !
  ! Open files
  !
  call dod_openfi(1_ip)
  !
  ! Initial variables
  !
  call dod_inivar(1_ip)
  !
  ! Read Dodeme data
  !
  call dod_readat()
  !
  ! Close data file
  !
  call dod_openfi(2_ip)
  !
  ! Write info
  !
  call dod_outinf()

end subroutine dod_turnon
