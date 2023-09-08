subroutine sld_outlat(itask)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_outlat
  ! NAME 
  !    sld_outlat
  ! DESCRIPTION
  !    This routine writes info on the heat equation in latex format
  ! USED BY
  !    sld_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solidz
  implicit none
  integer(ip),   intent(in) :: itask
  character(300)            :: equat
  character(20)             :: lvisc
  integer(ip)               :: ierhs
  character(500)            :: cgnup

  if(kfl_latex==1.and.kfl_paral<=0) then

     select case(itask)

     case(1)

     case(2)

     end select

  end if

end subroutine sld_outlat

