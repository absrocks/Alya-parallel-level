subroutine skyplu(neqns,lpdof,in_up,in_lo)

!-----------------------------------------------------------------------
!
! This routine obtains the positions IN_UP and IN_LO where the
! upper and lower matrices start.      
!
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip), intent(in)    :: neqns
  integer(ip), intent(in)    :: lpdof(neqns)
  integer(ip), intent(out)   :: in_up
  integer(ip), intent(inout) :: in_lo

  in_up = 1 + neqns
  in_lo = in_up + lpdof(neqns)

end subroutine skyplu
