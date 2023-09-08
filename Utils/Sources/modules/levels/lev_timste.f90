subroutine lev_timste
  !-----------------------------------------------------------------------
  !****f* Levels/lev_timste
  ! NAME 
  !    lev_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the level set 
  !    convection equation      
  ! USES
  !    lev_iniunk
  !    lev_updtss
  !    lev_updbcs
  !    lev_updunk
  !    lev_radvuf
  ! USED BY
  !    Levels
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_levels
  implicit none
  !
  ! Time step size 
  !
  if(kfl_stead_lev/=1) call lev_updtss()     

end subroutine lev_timste

