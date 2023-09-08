subroutine tur_begste()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_begste
  ! NAME 
  !    tur_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the turbulence
  !    equations
  ! USES
  !    tur_iniunk
  !    tur_updtss
  !    tur_updbcs
  !    tur_updunk
  !    tur_radvuf
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_turbul
  use def_kermod , only     :  kfl_prope
    implicit none
  !
  ! Variable density and viscosity
  !
  if (kfl_prope == 0_ip ) call tur_updpro()   ! With ker_proper this is not needed
  
  if(kfl_stead_tur/=1) then
     !
     ! Initial guess fo the turbulence variables: f(n,0,*) <-- f(n-1,*,*).
     !
     call tur_updunk(one)
     !
     ! Update boundary conditions
     !
     call tur_updbcs(TUR_BEFORE_TIME_STEP)

  end if
  !
  ! Calculate maximum mixing lenth for Apsey and Castro limitation (k eps)
  !
  if (kfl_lmaxi_tur==1.and.TUR_FAMILY_K_EPS) &
       call tur_maxlen()
end subroutine tur_begste
