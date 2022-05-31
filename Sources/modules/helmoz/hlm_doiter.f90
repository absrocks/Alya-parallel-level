subroutine hlm_doiter()

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_doiter.f90
  ! NAME 
  !    hlm_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the incompletely gauged
  !    coupled vector-scalar potential formulation of Maxwell's equations.
  !
  !    This formulation consists of two coupled equations that are solved 
  !    simultaneously: one vector equation which is modified vector Helm-
  !    holtz equation and one auxiliary scalar equation.
  ! USES
  !    hlm_begite
  !    hlm_solite
  !    hlm_endite
  ! USED BY
  !    Helmoz
  !-----------------------------------------------------------------------

  use def_parame
  use def_master
  use def_solver
  use def_helmoz

  implicit none

  !return

  !if(kfl_stead_hlm==0) then
  !call hlm_begite()
  !do while(kfl_goite_hlm==1)
   call hlm_solite()
  !   call hlm_endite(one)
  !end do
  !call hlm_endite(two)
  !end if

end subroutine hlm_doiter
