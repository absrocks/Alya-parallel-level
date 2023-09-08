subroutine chm_usrtem(T)
  !------------------------------------------------------------------------
  !****f* partis/chm_usrtem
  ! NAME 
  !    chm_usrtem
  ! DESCRIPTION
  !    This routine computes the temperature
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only  :  ip,rp
  use def_chemic, only  :  lawte_chm,temma_chm
  use def_master, only  :  cutim
  implicit none
  real(rp), intent(out) :: T

  if( lawte_chm == 0 ) then
     !
     ! Consant
     !
     T = temma_chm(1)

  else if( lawte_chm == 1 ) then
     !
     ! Linear funtion
     !
     T = temma_chm(1) + cutim * temma_chm(2)

  end if

end subroutine chm_usrtem
