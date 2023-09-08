subroutine lev_doiter
  !-----------------------------------------------------------------------
  !****f* Levels/lev_doiter
  ! NAME 
  !    lev_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the level set equation.
  ! USES
  !    lev_begite
  !    lev_solite
  !    lev_endite
  ! USED BY
  !    LEVELS
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_levels
  implicit none

  if(kfl_stead_lev==0) then
     call lev_begite
     do while(kfl_goite_lev==1)
        call lev_solite
        call lev_endite(one)
     end do
     call lev_endite(two)
  end if
  !
  ! Move mesh
  !
  !call lev_move_mesh_to_free_surface()

end subroutine lev_doiter
