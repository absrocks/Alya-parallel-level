subroutine nsi_updmsh(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_updmsh
  ! NAME
  !    updmsh
  ! DESCRIPTION
  !    This routine performs update of NST data. Depending on the 
  !    refinement strategy, some data is saved (to be used in the
  !    projections) and some arrays are deallocated (depending if 
  !    they change or not)
  ! USES
  ! USED BY
  !    nsi_newmsh
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use def_solver
  use mod_memchk
  implicit none
  integer(ip) :: itask

  select case(itask)

  case(one)

     if(kfl_algor_msh==1) then                  ! SHEAR_SLIP

     elseif(kfl_algor_msh==2) then              ! GLOBAL REFINEMENT


     end if

  case(two)

     if(kfl_algor_msh==1) then                  ! SHEAR_SLIP


     elseif(kfl_algor_msh==2) then              ! GLOBAL REFINEMENT

     end if

  end select

end subroutine nsi_updmsh
