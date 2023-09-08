subroutine nsa_begste
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_begste
  ! NAME 
  !    nsa_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step.  
  !    
  !    NASTAL module uses explicit schemes, until implicit ones
  !    were implemented. For that reason, nastal ALWAYS goes through 
  !    transients (kfl_timei_nsa > 0), even to solve pure stationary ones.
  !    
  !    The methods programmed are:
  !    
  !    Euler (for or back)             kfl_tisch_nsa=1
  !    BDF 2nd order                   kfl_tisch_nsa=2
  !    Crank-Nicolson                  kfl_tisch_nsa=3
  !    Improved Euler (Heun's)         kfl_tisch_nsa=30      r-k-3stages
  !    Crank-Nicholson                 kfl_tisch_nsa=31      r-k-3stages
  !    Kutta 3rd order                 kfl_tisch_nsa=41      r-k-3stages
  !    Kutta 4th order                 kfl_tisch_nsa=40      r-k-4stages
  !    
  ! USES
  !    nsa_iniunk
  !    nsa_updtss
  !    nsa_updbcs
  !    nsa_updunk
  ! USED BY
  !    Nastal
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain

  use def_nastal

  use mod_commdom_alya, only: INONE
  implicit none

  !call nsa_coupli(ITASK_BEGSTE)

  !
  ! Assign u(n,i-1,j) <-- u(n-1)
  !
  if(kfl_benme_nsa >= 200 .and. ittim > 1) then
     call nsa_apply_physics()
  end if 
  !
  ! Initial guess
  !
  call nsa_updunk(one)  ! u(,ITER_AUX) <-- u(,TIME_N) 

  !
  ! Update boundary conditions
  !
  call nsa_updbcs(one)

  call nsa_coupli(ITASK_BEGSTE)

end subroutine nsa_begste
