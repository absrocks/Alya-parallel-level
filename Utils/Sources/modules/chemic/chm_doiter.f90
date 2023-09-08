subroutine chm_doiter()
  !-----------------------------------------------------------------------
  !****f* partis/chm_doiter
  ! NAME 
  !    chm_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the transport equation.
  ! USES
  !    chm_begite
  !    chm_solite
  !    chm_endite
  ! USED BY
  !    Partis
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_chemic
  implicit none

  if(momod(modul) % kfl_stead==0) then
        call chm_begite()  
        do while(kfl_goite_chm==1 .and. kfl_reset /= 1)
           call chm_solite()                ! Solve PDE's
           call chm_endite(one)             ! Update unknowns
        end do
        call chm_endite(two)
  end if
end subroutine chm_doiter
