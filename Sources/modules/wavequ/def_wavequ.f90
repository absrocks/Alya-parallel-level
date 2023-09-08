module def_wavequ
  !-----------------------------------------------------------------------
  !    
  ! Heading for the wave equation routines
  !
  !-----------------------------------------------------------------------
  use def_kintyp
  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------
       
  real(rp),      parameter       :: &
       zewav = epsilon(1.0_rp)
  integer(ip),   parameter       :: &
       ncoef_wav=10,                & ! # coefficient for properties
       nsour_wav=10                   ! # coefficient for source term

  !------------------------------------------------------------------------
  ! Physical problem
  !------------------------------------------------------------------------
!--BEGIN REA GROUP
  integer(ip)                    :: &
       kfl_sourc_wav                  ! If there is a source term
  real(rp)                       :: &
       sourc_wav(nsour_wav)           ! Parameter for the source term
  !
  ! Physical properties
  !
  integer(ip)                    :: &
       nmate_wav                      ! # of materials
  integer(ip),  allocatable ::      &
       lawde_wav(:),                & ! Law for rho
       lawka_wav(:)                   ! Law for k
  integer(ip),  pointer          :: &
       lmate_wav(:)                   ! Medium material (element-wise)
  real(rp),     allocatable      :: &
       densi_wav(:,:),              & ! Density (rho)
       kappa_wav(:,:)                 ! Compressibility kappa (k)

  !------------------------------------------------------------------------
  ! Numerical problem
  !------------------------------------------------------------------------
  integer(ip)                    :: &
       kfl_tiacc_wav,               & ! Temporal accuracy
       kfl_tisch_wav,               & ! Time integration scheme
       kfl_normc_wav,               & ! Norm of convergence
       kfl_timet_wav,               & ! Time treatment
       kfl_subgs_wav,               & ! Subgrid scale
       kfl_massm_wav,               & ! Type of mass matrix
       neule_wav,                   & ! Number of Euler time steps 
       miinn_wav                      ! Max # of iterations
  real(rp)                       :: &
       safet_wav,                   & ! Safety factor for time step
       sstol_wav,                   & ! Steady state tolerance
       cotol_wav,                   & ! Convergence tolerance
       nebet_wav,                   & ! Beta of Newmark scheme
       negam_wav                      ! Gamma of Newmark scheme
 
  !------------------------------------------------------------------------
  ! Output and Postprocess
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Boundary conditions
  !------------------------------------------------------------------------
  integer(ip)                    :: &
       kfl_onnod_wav,               & ! If there exists bc on boundaries
       kfl_onbou_wav,               & ! If there exists bc on boundaries
       kfl_absor_wav                  ! If there exists absorbing bc
  integer(ip), pointer           :: &
       kfl_fixno_wav(:),            & ! Nodal fixity 
       kfl_fixbo_wav(:)               ! Boundary fixity 
  real(rp),    pointer           :: &
       bvess_wav(:)                   ! Essential bc values
!--END REA GROUP
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------
  integer(ip)                    :: &
       ncomp_wav,                   & ! Number of components of the temperature
       kfl_tiaor_wav,               & ! Original time accuracy
       kfl_stead_wav,               & ! Steady-state has been reached 
       kfl_goite_wav,               & ! Keep iterating
       nzmat_wav,                   & ! Matrix size
       nzrhs_wav                      ! RHS size
  real(rp)                       :: &
       dtcri_wav,                   & ! Critical time step
       cpuit_wav,                   & ! CPU time per iteration
       dtinv_wav,                   & ! 1/dt
       pabdf_wav(10),               & ! BDF factors
       wamin_wav,                   & ! Minimum wave amplitude
       wamax_wav,                   & ! Maximum wave amplitude
       resid_wav                      ! Residual for outer iterations
  real(rp),    pointer           :: &
       vmass_wav(:),                & ! Modified matrix when using absorbing bc
       wavve_wav(:,:),              & ! Wave velocity
       wavac_wav(:,:)                 ! Wave acceleration
  type(r2p), pointer             :: &
       wasgs(:)                      ! Wave amplitud subgrid scale  

end module def_wavequ
