subroutine chm_begite()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_begite
  ! NAME 
  !    chm_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the transport
  !    equation
  ! USES
  !    chm_tittim
  !    chm_updbcs
  !    chm_inisol
  !    chm_updunk
  ! USED BY
  !    chm_doiter
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_chemic
  use def_solver
  use def_kermod
  use mod_ker_proper 
  use mod_messages, only : livinf
  implicit none
  !
  ! Initializations
  !
  kfl_goite_chm = 1 
  itinn(modul)  = 0
  kfl_under_chm = 0
  kfl_overs_chm = 0

  if(itcou==1 .and. kfl_model_chm /=4) call chm_tistep()
  call livinf(15_ip,' ',modul)

  if (kfl_model_chm /=4) then
     !
     ! Update boundary conditions
     !
     call chm_updbcs(ITASK_BEGITE)

     !
     ! Obtain the initial guess for inner iterations
     !
     call chm_updunk(ITASK_BEGITE)
  end if
end subroutine chm_begite
    
