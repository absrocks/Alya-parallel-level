subroutine shape0(nnode,shapf,ierro)

  !-----------------------------------------------------------------------
  !
  ! This routine evaluates shape functions and their first derivates 
  ! for 0-d continuous with 1 nodes
  !
  !-----------------------------------------------------------------------

  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: nnode
  integer(ip), intent(out) :: ierro
  real(rp),    intent(out) :: shapf(nnode)

  ierro    = 0 
  shapf(1) = 1.0_rp

end subroutine shape0
