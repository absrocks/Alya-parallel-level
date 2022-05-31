subroutine shafal(&
     posgp,pdime,pnode,pgaus,ntens,gpsha,gpder,gphes,ierro)
  !-----------------------------------------------------------------------
  !
  !    This routine evaluates shape functions and their gpderatives
  !    for linear and quadratic isoparametric elements
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)    :: pdime
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: ntens
  integer(ip), intent(inout) :: ierro
  real(rp),    intent(in)    :: posgp(pdime,pgaus)
  real(rp),    intent(out)   :: gpsha(pnode,pgaus)
  real(rp),    intent(out)   :: gpder(max(1_ip,pdime),pnode,pgaus)
  real(rp),    intent(out)   :: gphes(max(1_ip,ntens),pnode,pgaus)
  integer(ip)                :: igaus

  ierro = 0
  gpsha = 0.0_rp
  gpder = 0.0_rp
  gphes = 0.0_rp

  do igaus = 1,pgaus

     if( pdime == 0 ) then
        call shape0(pnode,gpsha(1,igaus),ierro)

     else if( pdime == 1 ) then
        call shape1(posgp(1,igaus),pnode,gpsha(1,igaus),gpder(1,1,igaus),&
             gphes(1,1,igaus),ierro)

     else if( pdime == 2 ) then
        call shape2(posgp(1,igaus),posgp(2,igaus),pnode,gpsha(1,igaus),gpder(1,1,igaus),&
             gphes(1,1,igaus),ierro)

     else if( pdime == 3 ) then
        call shape3(posgp(1,igaus),posgp(2,igaus),posgp(3,igaus),pnode,gpsha(1,igaus),&
             gpder(1,1,igaus),gphes(1,1,igaus),ierro)

     end if

  end do

end subroutine shafal
