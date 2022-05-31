subroutine tem_begste()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_begste
  ! NAME 
  !    tem_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step 
  ! USES
  !    tem_iniunk
  !    tem_updtss
  !    tem_updbcs
  !    tem_updunk
  !    tem_radvuf
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  implicit none
  
  if(kfl_stead_tem/=1) then     
     !
     ! Initial guess fo the temperature: T(n,0,*) <-- T(n-1,*,*).
     !
     call tem_updunk(ITASK_BEGSTE)
     !
     ! Update boundary conditions
     !
     call tem_updbcs(ITASK_BEGSTE)
     !
     ! Compute view factors if radiation is considered
     !
     if(kfl_radia_tem==1.and.ittim==1) call tem_radvuf()
     !
     ! Coupling with dynamic solver
     !
     call tem_dyncou(1_ip)   

  end if

  call tem_coupli(ITASK_BEGSTE)
    
end subroutine tem_begste

