subroutine chm_doiter()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_doiter
  ! NAME 
  !    chm_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the transport equation.
  ! USES
  !    chm_begite
  !    chm_solite
  !    chm_endite
  ! USED BY
  !    Chemic
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_chemic
  use mod_chm_explicit,    only : chm_explicit_solution
  use mod_chm_rk_explicit, only : chm_rk_explicit_solution
  use mod_chm_levSet,      only : chm_pseudo_time_reinit_levSet
  implicit none

  if(momod(modul)%kfl_stead == 0) then

        call chm_begite()  

        do while( kfl_goite_chm==1 )

           if (kfl_tisch_chm == 3) then
              call chm_explicit_solution()         ! Solve AB scheme
           else
              if(kfl_spray_chm == 2) then
                 call chm_pseudo_time_reinit_levSet()

                 call chm_rk_explicit_solution()   ! Solve surface density PDE

              else
                 call chm_rk_explicit_solution()   ! Solve RK scheme
              end if
           endif

        end do

        call chm_endite(ITASK_ENDITE)

  end if


end subroutine chm_doiter
