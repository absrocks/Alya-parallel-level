subroutine chm_begste()
  !-----------------------------------------------------------------------
  !****f* partis/chm_begste
  ! NAME 
  !    chm_begste
  ! DESCRIPTION
  !    This routine prepares a new time step
  ! USES
  !    chm_updunk
  ! USED BY
  !    partis
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_memory
  use mod_chm_finiteRate, only : chm_getProp_finiteRate
  implicit none

  if (kfl_model_chm /= 4) then
     if(  momod(modul) % kfl_stead /=1) then     
        !
        ! Initial guess: c(n,0,*) <-- c(n-1,*,*).
        !
        call chm_updunk(ITASK_BEGSTE)
     end if
  end if

  if (kfl_model_chm == 1) then

     !
     ! Read flamelet table for gas phase
     !
     if (kfl_spray_chm == 0 .or. ( kfl_spray_chm /= 0 .and. kfl_premix_chm == 0)) then
        if (kfl_ufpv_chm > 0) then
           call chm_post_scalar_dissipation_rate(47_ip)
        endif

        if (kfl_lookg_chm > 0) then
            call chm_gp_reatab()
        else
            call chm_reatab()
        endif
     end if
 
  elseif (kfl_model_chm == 3) then
       !
       ! Calculate transport properties
       !
       call chm_getProp_finiteRate()

  endif
  
end subroutine chm_begste

