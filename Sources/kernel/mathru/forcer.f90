subroutine forcer()

  !------------------------------------------------------------------------
  !
  ! This function results in an error so that alya stops and the traceback tells where you are
  !
  !------------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp
  implicit none

  real(rp)                 :: hhhhh

  call cputim(hhhhh)

  hhhhh = sin(hhhhh)
  hhhhh = sqrt(hhhhh-5.0)
  
end subroutine forcer
 
