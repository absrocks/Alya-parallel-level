subroutine elsest_boxcoo(ndime,box_nr,box_coord)
  ! returns the box coordinates in i j k of a given box number
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip),intent(in)          :: ndime,box_nr
  integer(ip),intent(out)         :: box_coord(3)
  real(rp)                        :: temp,temp2,temp3,box_nr_real
  real(rp)                        :: box_coord_real(3),nboxxreal(3)

  box_coord_real = real(box_coord,rp)
  box_nr_real    = real(box_nr,rp)
  nboxxreal      = real(nboxx,rp)

  if(ndime==3) then
     temp              = box_nr_real/(nboxxreal(1)*nboxxreal(2)+0.00001_rp)
     box_coord(3)      = int(temp,ip)+1
     box_coord_real(3) = real(box_coord(3),rp)
     temp2             = box_nr_real - int(temp,ip)*nboxxreal(1)*nboxxreal(2)
     temp3             = temp2/(nboxxreal(1)+0.00001_rp)
     box_coord(2)      = int(temp3,ip)+1
     box_coord_real(2) = real(box_coord(2),rp)

     box_coord(1)      =   int((box_nr_real-( (box_coord_real(3)- 1.0_rp)*nboxxreal(1)*nboxxreal(2))) &
          &              - ((box_coord_real(2)-1.0_rp) * nboxxreal(1)),ip)
  else
     temp              = box_nr_real/(nboxxreal(1)+0.00001_rp)
     box_coord(2)      = int(temp,ip) + 1
     box_coord_real(2) = real(box_coord(2),rp)
     temp2             = box_nr_real - int(temp,ip)*nboxxreal(1)
     box_coord(1)      = int(temp2,ip)
  end if

end subroutine elsest_boxcoo
