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
  implicit none

  if(  momod(modul) % kfl_stead /=1) then     
     !
     ! Initial guess: c(n,0,*) <-- c(n-1,*,*).
     !
     call chm_updunk(one)
     call chm_upwmea(2_ip) ! Mean molecular weight
  end if
  !
  ! Meteo model
  !
  if( kfl_model_chm == 2 ) call chm_begmet()

end subroutine chm_begste

