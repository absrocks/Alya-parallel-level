subroutine elsest_neighb(&
     ndime,nboxx,curr_box_coor,advanzable_less,advanzable_more,ilevel)
  !
  ! Look for neighboring boxes
  !
  use def_elsest, only      : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: nboxx(ndime),curr_box_coor(ndime)
  integer(ip), intent(in)  :: ilevel
  integer(ip), intent(out) :: advanzable_less(ndime),advanzable_more(ndime)
  integer(ip)              :: ii

  do ii=1,ndime
     if(curr_box_coor(ii)<=ilevel) then
        advanzable_less(ii) = -1
     else
        advanzable_less(ii) =  1
     end if
     if(curr_box_coor(ii)>nboxx(ii)-ilevel) then
        advanzable_more(ii) = -1
     else
        advanzable_more(ii) =  1
     end if
  end do

end subroutine elsest_neighb
