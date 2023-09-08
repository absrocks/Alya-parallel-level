subroutine chm_turnof()
  !-----------------------------------------------------------------------
  !****f* partis/chm_turnof
  ! NAME 
  !    chm_turnof
  ! DESCRIPTION
  !    This routine closes the run for the temperature equation
  ! USES
  !    chm_outcpu
  ! USED BY
  !    Partis
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_chemic
  implicit none
  !
  ! Results
  !
  call chm_sizedi()
  call chm_sizusr()
  !
  ! Close files
  !
  call chm_openfi(4_ip)
  !
  ! Output CPU times
  !
  call chm_outcpu()

  !-------------------------------------------------------------------
  !
  ! Interpolation from coarse to fine mesh
  !
  !-------------------------------------------------------------------
  !
  if(kfl_meshi_chm /= 0_ip) call chm_coarfine(2_ip)

end subroutine chm_turnof

