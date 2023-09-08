subroutine got_timste
  !-----------------------------------------------------------------------
  !****f* Gotita/got_timste
  ! NAME 
  !    got_timste
  ! DESCRIPTION
  !    This routine computes the time step
  ! USES
  !    got_iniunk
  !    got_updtss
  !    got_updbcs
  !    got_updunk
  ! USED BY
  !    Gotita
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  implicit none
  !
  ! Latex output format
  !
  if(ittim==0) call got_outlat(one)
  !
  ! For the first time step, set up the initial condition and initialize
  !
  if(ittim==0) then
     kfl_stead_got = 0
  end if
  !
  ! Time step size
  !
  if(kfl_stead_got/=1) call got_updtss()      
  
end subroutine got_timste
