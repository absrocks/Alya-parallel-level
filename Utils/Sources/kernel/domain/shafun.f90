subroutine shafun(&
     posgp,ndime,nnode,ntens,shapf,deriv,heslo,ierro)

  !-----------------------------------------------------------------------
  !
  !    This routine evaluates shapf functions and their derivatives
  !    for linear and quadratic isoparametric elements
  !
  !-----------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  implicit none

  integer(ip), intent(in)    :: ndime,nnode,ntens
  integer(ip), intent(inout) :: ierro
  real(rp),    intent(in)    :: posgp(ndime)
  real(rp),    intent(out)   :: shapf(nnode),deriv(max(1_ip,ndime),nnode)
  real(rp),    intent(out)   :: heslo(max(1_ip,ntens),nnode)
  !
  ! Initializations
  !
  shapf = 0.0_rp
  deriv = 0.0_rp
  ierro = 0
  if( ndime /= 0 ) heslo = 0.0_rp
  !
  ! Evaluation of the shapf functions
  !
  if( ndime == 0 ) then
     call shape0(nnode,shapf,ierro)
  else if( ndime == 1 ) then
     call shape1(posgp(1),nnode,shapf,deriv,heslo,ierro)
  else if( ndime == 2 ) then
     call shape2(posgp(1),posgp(2),nnode,shapf,deriv,heslo,ierro)
  else if( ndime == 3 ) then
     call shape3(posgp(1),posgp(2),posgp(3),nnode,shapf,deriv,heslo,ierro)
  end if

end subroutine shafun
