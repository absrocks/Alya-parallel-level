subroutine chm_newmsh()
  !-----------------------------------------------------------------------
  !****f* partis/chm_newmsh
  ! NAME 
  !    chm_newmsh
  ! DESCRIPTION
  !    This routine projects data from the old to the new mesh
  ! USES
  ! USED BY
  !    Partis
  !***
  !-----------------------------------------------------------------------
  use def_master
  implicit none

  if(kfl_algor_msh==1) then                  ! SHEAR_SLIP

  else if(kfl_algor_msh==2) then              ! GLOBAL REFINEMENT

  end if

end subroutine chm_newmsh
