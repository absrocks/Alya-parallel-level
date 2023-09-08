program test_elmgeo_natural_coordinates_2d
  
  use def_kintyp, only : ip,rp,lg
  use mod_elmgeo, only : elmgeo_natural_coordinates
  use mod_elmgeo, only : element_type
  use def_elmtyp
  use unitt
  implicit none
  integer(ip), parameter :: mnode = 64
  integer(ip), parameter :: ndime = 3

  integer(ip)            :: pnode
  integer(ip)            :: pelty
  real(rp)               :: elcod(ndime,mnode)
  real(rp)               :: coloc(3)
  
  real(rp)               :: coglo(3)
  real(rp)               :: shapf(mnode)
  real(rp)               :: deriv(ndime,mnode)
  real(rp)               :: toler
  real(rp)               :: xx(3),xmax,error
  integer(ip)            :: inode,ifoun
  !
  ! Select my test
  !
  call my_test(&
       2_ip,pelty,pnode,elcod,coglo,toler)
  !
  ! Inclusion test
  !
  call elmgeo_natural_coordinates(             &
       ndime,pelty,pnode,elcod,shapf,deriv,    &
       coglo,coloc,ifoun,toler)
  !
  ! Test coordinates are ok 
  !
  xx = 0.0_rp
  do inode = 1,pnode
     xx(1:ndime) = xx(1:ndime) + shapf(inode) * elcod(1:ndime,inode)
  end do
  xmax  = sqrt( maxval(coglo)-minval(coglo) )
  error = sqrt( dot_product(xx(1:ndime)-coglo(1:ndime),xx(1:ndime)-coglo(1:ndime)) )

  !call assert_true(error/xmax>1.0e-12_rp,'ELMGEO_NATURAL_COORDINATES','1')
  call assert_true(.false.,'ELMGEO_NATURAL_COORDINATES','1')
  
contains

  subroutine my_test(itest,pelty,pnode,elcod,coglo,toler1)
    
    integer(ip), intent(in)  :: itest
    integer(ip), intent(out) :: pelty
    integer(ip), intent(out) :: pnode
    real(rp),    intent(out) :: elcod(ndime,mnode)
    real(rp),    intent(out) :: coglo(3)
    real(rp),    intent(out) :: toler1

    toler1 = 1.0e-6_rp
    
    select case ( itest )

    case ( 1_ip )
       !
       ! Simple QUA04
       !
       pelty      = QUA04
       coglo(1)   = 0.5_rp
       coglo(2)   = 0.5_rp
       elcod(1,1) = 0.0_rp
       elcod(2,1) = 0.0_rp
       elcod(1,2) = 1.0_rp
       elcod(2,2) = 0.0_rp
       elcod(1,3) = 1.0_rp
       elcod(2,3) = 1.0_rp
       elcod(1,4) = 0.0_rp
       elcod(2,4) = 1.0_rp

    case ( 2_ip )
       !
       ! PYR05
       !
       pelty        = PYR05
       coglo        = (/ 9.874581150000000E-002_rp , 0.269321775900000_rp , -0.301004499000000_rp /)     
       elcod(1:3,1) = (/ 9.888322240000000E-002_rp , 0.269127140100000_rp , -0.301096436900000_rp /)
       elcod(1:3,2) = (/ 9.902521950000000E-002_rp , 0.269250333900000_rp , -0.300936568800000_rp /)     
       elcod(1:3,3) = (/ 9.903627550000001E-002_rp , 0.269340224000000_rp , -0.300967480300000_rp /)     
       elcod(1:3,4) = (/ 9.892776609999999E-002_rp , 0.269234288200000_rp , -0.301114527300000_rp /)     
       elcod(1:3,5) = (/ 9.874581150000000E-002_rp , 0.269321775900000_rp , -0.301004499000000_rp /)    

    end select

    pnode = element_type(pelty) % number_nodes
    
  end subroutine my_test
    
end program test_elmgeo_natural_coordinates_2d

 
