
subroutine tem_doiter()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_doiter
  ! NAME 
  !    tem_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the temperature equation.
  ! USES
  !    tem_begite
  !    tem_solite
  !    tem_endite
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_temper
  implicit none

  if( kfl_stead_tem == 0 ) then
     call tem_begite()
     do while( kfl_goite_tem == 1 )
        call tem_solite()
        call tem_endite(one)
     end do
     call tem_endite(two)
  end if

end subroutine tem_doiter
