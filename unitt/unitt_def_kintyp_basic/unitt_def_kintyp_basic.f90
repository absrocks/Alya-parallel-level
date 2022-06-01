program unitt_def_kintyp_basic 
  
  use def_kintyp, only : ip, rp

  implicit none

  if (ip == 4 .or. ip == 8) then
    print *, "ip value is ", ip
  else
    print *, "ip value is invalid: ", ip
    stop 1
  end if
  if (rp == 4 .or. rp == 8 .or. rp == 16) then
    print *, "rp value is ", ip
  else
    print *, "rp value is invalid: ", rp
    stop 1
  end if

end program unitt_def_kintyp_basic
