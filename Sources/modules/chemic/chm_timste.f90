subroutine chm_timste()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_timste
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
  !    Chemic
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use def_solver
  implicit none

  !
  ! Time step size 
  !
  if(momod(modul) % kfl_stead/=1) call chm_updtss()     

end subroutine chm_timste
