subroutine rupoin(ngaus,weigp,ierro)
  !-----------------------------------------------------------------------
  !****f* Domain/rupoin
  ! NAME
  !    rupoin
  ! DESCRIPTION
  !     This routine sets up the integration constants
  ! OUTPUT
  !    
  ! USED BY
  !    rulepw
  !***
  !-----------------------------------------------------------------------
  use  def_kintyp, only    :  ip,rp
  implicit none
  integer(ip), intent(in)  :: ngaus
  integer(ip), intent(out) :: ierro
  real(rp),    intent(out) :: weigp(ngaus)

  ierro    = 0
  weigp(1) = 1.0_rp

end subroutine rupoin
