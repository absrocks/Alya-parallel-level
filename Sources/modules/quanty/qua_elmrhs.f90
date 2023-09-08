subroutine qua_elmrhs(pgaus,gpmas,gprhs)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_elmrhs
  ! NAME
  !   qua_elmrhs
  ! DESCRIPTION
  !    Compute right-hand side of the equation
  ! OUTPUT 
  !    GPRHS
  ! USES
  ! USED BY
  !    qua_elmope
  !***
  !-----------------------------------------------------------------------
  use def_quanty

  implicit none
  integer(ip), intent(in)    :: pgaus
  real(rp),    intent(in)    :: gpmas(pgaus)
  real(rp),    intent(out)   :: gprhs(pgaus)
  integer(ip)                :: igaus

  !
  ! Source term: Q
  !
  do igaus=1,pgaus
     gprhs(igaus)=gpmas(igaus)
  end do

end subroutine qua_elmrhs
