subroutine elsest_boxnum(ndime,nboxx,box_coord,box_nr)
  !
  ! Returns the box number box_nr of a given box with coordinates in i j k
  !
  use def_elsest, only      : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: nboxx(ndime),box_coord(ndime)
  integer(ip), intent(out) :: box_nr

  if(ndime==2) then
     box_nr=(box_coord(2)-1)*nboxx(1) + box_coord(1)
  else
     box_nr=(box_coord(3)-1)*(nboxx(1)*nboxx(2)) + (box_coord(2)-1)*nboxx(1) + box_coord(1)
  end if

end subroutine elsest_boxnum
