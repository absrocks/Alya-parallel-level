subroutine rulepw(ndime,ngaus,nrule,posgp,weigp,ierro)
  !-----------------------------------------------------------------------
  !****f* domain/rulepw
  ! NAME
  !    rulepw
  ! DESCRIPTION
  !    This routine sets up the integration constants
  ! OUTPUT
  !    An n-point Gaussian quadrature rule is a quadrature rule constructed 
  !    to yield an exact result for polynomials of degree 2n-1
  !    POSGP(NDIME,NGAUS) ... Local coordinates of Gauss points
  !    WEIGP(NGAUS) ......... Weights of Gauss points
  ! USES
  !    ruqope
  !    ruqclo
  !    rutope
  !    rutclo
  !    rupope
  !    rupclo
  !    ruyope
  !    ruyclo
  ! USED BY
  !    domain
  !***
  !----------------------------------------------------------------------- 
  use      def_kintyp
  implicit none
  integer(ip), intent(in)    ::  ndime
  integer(ip), intent(in)    ::  ngaus
  integer(ip), intent(in)    ::  nrule
  integer(ip), intent(inout) ::  ierro
  real(rp),    intent(out)   ::  posgp(max(1_ip,ndime),ngaus), weigp(ngaus)

  ierro = 0
  posgp = 0.0_rp
  weigp = 0.0_rp

  select case(nrule)

  case(-2_ip) 
     call rupoin(ngaus,weigp,ierro)             ! point
  case(-1_ip) 
     call rubope(ngaus,posgp,weigp,ierro)       ! Line:         open
  case( 0_ip) 
     call rubclo(ngaus,posgp,weigp,ierro)       ! Line:         closed
  case( 1_ip) 
     call ruqope(ndime,ngaus,posgp,weigp,ierro) ! Quad/Hexa:    open
  case( 2_ip) 
     call ruqclo(ndime,ngaus,posgp,weigp,ierro) ! Quad/Hexa:    closed
  case( 3_ip) 
     call rutope(ndime,ngaus,posgp,weigp,ierro) ! Tria/Tetra:   open 
  case( 4_ip) 
     call rutclo(ndime,ngaus,posgp,weigp,ierro) ! Tria/Tetra:   closed 
  case( 5_ip) 
     call rupope(ndime,ngaus,posgp,weigp,ierro) ! Prism(Penta): open 
  case( 6_ip) 
     call rupclo(ndime,ngaus,posgp,weigp,ierro) ! Prism(Penta): closed
  case( 7_ip) 
     call ruyope(ndime,ngaus,posgp,weigp,ierro) ! Pyramid:      open
  case( 8_ip) 
     call ruyclo(ndime,ngaus,posgp,weigp,ierro) ! Pyramid:      closed

  end select

end subroutine rulepw
