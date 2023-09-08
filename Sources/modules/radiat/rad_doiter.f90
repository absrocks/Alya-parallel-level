subroutine rad_doiter()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_doiter
  ! NAME 
  !    rad_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the radiation equation.
  ! USES
  !    rad_begite
  !    rad_solite
  !    rad_endite
  ! USED BY
  !    Radiat
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_radiat
  implicit none

  call rad_begite()
  do while(kfl_goite_rad==1)
     call rad_solite()
     call rad_endite(one)
  end do
  call rad_endite(two)

end subroutine rad_doiter
