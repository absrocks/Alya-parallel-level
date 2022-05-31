subroutine chm_turnof()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_turnof
  ! NAME 
  !    chm_turnof
  ! DESCRIPTION
  !    This routine closes the run for the chemic module
  ! USES
  !    chm_outcpu
  ! USED BY
  !    Chemic
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_chemic
  implicit none

  !
  ! Close files
  !
  call chm_openfi(4_ip)

  !
  ! Output CPU times
  !
  call chm_outcpu()

end subroutine chm_turnof

