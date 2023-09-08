subroutine elsest_nextbo(ndime,ithre, &
     curr_box_coor,advanzable_less, &
     & advanzable_more,ilevel,chbox_nr)

  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), intent(in)                      :: ndime,ithre
  integer(ip), intent(in)                      :: ilevel
  integer(ip), intent(in)                      :: advanzable_less(3),advanzable_more(3)
  integer(ip), intent(in)                      :: chbox_nr
  integer(ip), intent(in)                      :: curr_box_coor(3)
  integer(ip)                                  :: temp_box_coor(3)
  integer(ip)                                  :: counter_min(3),counter_max(3)
  integer(ip)                                  :: i,j,k,counter,temp_box_nr
  integer(ip)                                  :: i_t,j_t, k_t,m
  logical                                      :: box_checked

  nexbo(:,ithre)=0
  counter_min=0
  counter_max=0
  do i=1,ndime
     if (advanzable_less(i) .eq. 1) then
        counter_min(i) = -ilevel
     end if
     if (advanzable_more(i) .eq. 1) then
        counter_max(i) = ilevel
     end if
  end do

  counter=0
  if(ndime==3) then
     do i=counter_min(3),counter_max(3)
        do j=counter_min(2),counter_max(2)
           do k=counter_min(1),counter_max(1)
              box_checked = .false.
              temp_box_coor(1)=curr_box_coor(1)+k
              temp_box_coor(2)=curr_box_coor(2)+j
              temp_box_coor(3)=curr_box_coor(3)+i
              call elsest_boxnum(ndime,nboxx,temp_box_coor,temp_box_nr)
              do m=1,chbox_nr
                 if(chbox(m,ithre) .eq. temp_box_nr) then
                    box_checked = .true.
                 end if
              end do
              if(box_checked.eqv..false.) then
                 counter=counter+1
                 nexbo(counter,ithre) = temp_box_nr
                 if (i /= 0) then
                    i_t = i/abs(i)
                 else
                    i_t = i
                 end if
                 if (j /= 0) then
                    j_t = j/abs(j)
                 else
                    j_t = j
                 end if
                 if (k /= 0) then
                    k_t = k/abs(k)
                 else
                    k_t = k
                 end if
              end if
           end do
        end do
     end do
  else if(ndime==2) then
     do j=counter_min(2),counter_max(2)
        do k=counter_min(1),counter_max(1)
           box_checked = .false.
           temp_box_coor(1)=curr_box_coor(1)+k
           temp_box_coor(2)=curr_box_coor(2)+j
           call elsest_boxnum(ndime,nboxx,temp_box_coor,temp_box_nr)
           do m=1,chbox_nr
              if(chbox(m,ithre) .eq. temp_box_nr) then
                 box_checked = .true.
              end if
           end do
           if(box_checked.eqv..false.) then
              counter=counter+1
              nexbo(counter,ithre) = temp_box_nr
              if (i /= 0) then
                 i_t = i/abs(i)
              else
                 i_t = i
              end if
              if (j /= 0) then
                 j_t = j/abs(j)
              else
                 j_t = j
              end if
              if (k /= 0) then
                 k_t = k/abs(k)
              else
                 k_t = k
              end if
           end if
        end do
     end do
  end if
end subroutine elsest_nextbo
