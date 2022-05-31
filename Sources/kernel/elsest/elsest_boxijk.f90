subroutine elsest_boxijk(ndime,nboxx,point_x,curr_box_coor,comin,comax)
  !
  ! Look for the box curr_box_coor)(i,j,k) containing point point_x
  !
  use def_elsest, only      : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: nboxx(ndime)
  real(rp),    intent(in)  :: point_x(ndime),comin(ndime),comax(ndime)
  integer(ip), intent(out) :: curr_box_coor(ndime)
  integer(ip)              :: idime

  do idime = 1,ndime
     curr_box_coor(idime) = int( ( (point_x(idime) - comin(idime)) / (comax(idime) - comin(idime)) ) * real(nboxx(idime),rp) , ip ) + 1    
     curr_box_coor(idime) = min( curr_box_coor(idime) , nboxx(idime) )
  end do

end subroutine elsest_boxijk
