subroutine chm_timste()
  !-----------------------------------------------------------------------
  !****f* partis/chm_timste
  ! NAME 
  !    chm_begste
  ! DESCRIPTION
  !    This routine computes the new time step
  ! USES
  !    chm_iniunk
  !    chm_updtss
  !    chm_updbcs
  !    chm_updunk
  !    chm_radvuf
  ! USED BY
  !    Partis
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use def_solver
  implicit none
  !
  ! Update fields
  !
  call chm_updfie()
  if (kfl_model_chm==4) then
     call chm_omegak(1_ip,nspec_chm)     ! Mass source rates
  endif
  !
  ! Time step size 
  !
  if(momod(modul) % kfl_stead/=1) call chm_updtss()     

end subroutine chm_timste
