
subroutine Endopt()
  !-----------------------------------------------------------------------
  !****f* master/Endopt
  ! NAME
  !    Endopt
  ! DESCRIPTION
  !    This routine accepts/rejects the update candidate for design vars, 
  !    updates the step length in a line-search and ends the optimization
  !    if convergence is achieved.
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  implicit none

  if(kfl_servi(ID_OPTSOL)==0) then
     kfl_goopt=0
  end if

  call Optsol(ITASK_ENDOPT)

end subroutine Endopt
