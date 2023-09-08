subroutine got_outlat(itask)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_outlat
  ! NAME 
  !    got_outlat
  ! DESCRIPTION
  !    This routine writes info on the incomcdropible Navier-Stokes 
  !    equations in latex format
  ! USED BY
  !    got_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_gotita
  implicit none
  integer(ip), intent(in) :: itask

  if(kfl_latex==1.and.kfl_paral<=0) then


  end if

end subroutine got_outlat
      
