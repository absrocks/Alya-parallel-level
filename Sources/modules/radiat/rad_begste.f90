subroutine rad_begste()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_begste
  ! NAME 
  !    rad_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step 
  ! USES
  !    rad_iniunk
  !    rad_updtss
  !    rad_updbcs
  !    rad_updunk
  !    rad_radvuf
  ! USED BY
  !    Radiat
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  implicit none
  integer(ip), save :: ipass=0

  if(ipass==0) then
     !
     ! First time we pass here: do not use ittim as this can be 
     ! a restart run 
     !
     ipass=1

  end if
  !
  ! Initial guess fo the radiation: G(n,0,*) <-- G(n-1,*,*).
  !
  call rad_updunk(one)
  !
  ! Update boundary conditions
  !
  call rad_updbcs(one)   !!F Need to put here some of the update of the BCS from previous iterations...
  !
  ! Coupling with dynamic solver  
  !
  !call rad_dyncou(1_ip)
  
end subroutine rad_begste

