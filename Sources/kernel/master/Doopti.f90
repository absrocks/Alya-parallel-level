
subroutine Doopti()
  !-----------------------------------------------------------------------
  !****f* master/Doopti
  ! NAME
  !    Doopti
  ! DESCRIPTION
  !    This routine calculates a descent direction and an update candidate
  !    for the design variables.
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  implicit none
  !
  ! End a time step for each module
  !
  call moduls(ITASK_DOOPTI)

  !if(kfl_servi(ID_OPTSOL)==1) then
     call Optsol(ITASK_DOOPTI)
  !end if

  end subroutine Doopti
