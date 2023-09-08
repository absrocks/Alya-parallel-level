subroutine got_begste
  !-----------------------------------------------------------------------
  !****f* Gotita/got_begste
  ! NAME 
  !    got_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the incomcdropible NS
  !    equations.
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
  integer(ip), save :: ipass=0

  if(ipass==0) then
     !
     ! First time we pass here: do not use ittim as this can be 
     ! a restart run 
     !
     ipass=1
     call got_inivar(two)
  end if

  if(kfl_stead_got/=1)  then     
     !
     ! Initial guess fo the droplet velocity: u(n,0,*) <-- u(n-1,*,*).
     !
     call got_updunk(one)
     !
     ! Initialize some variables before starting the time step
     !
     call got_inivar(three)
     !
     ! Stop artificial viscosity and shock capturing
     !
     if(kfl_artif_got/=0) then
        if(ittim>=itart_got) kfl_artif_got=0
     end if
     if(kfl_shocm_got/=0) then
        if(ittim>=itshm_got) kfl_shocm_got=0
     end if
     if(kfl_shocc_got/=0) then
        if(ittim>=itshc_got) kfl_shocc_got=0
     end if

     if(ittim>=100) then
        !solve(1)%kfl_algso=-1
     end if
  end if

end subroutine got_begste
