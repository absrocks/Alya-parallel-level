module def_gotita
  !------------------------------------------------------------------------
  !    
  ! Heading for the incomcdropible GOTITA routines
  !
  !------------------------------------------------------------------------
  use def_kintyp
  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------
  integer(ip),   parameter :: &
       lun_bound_got = 110, lun_stasg_got = 112, lun_cvgsg_got = 113
  character(150)           :: &
       fil_conve_got                 ! Convergence file name
  integer(ip),   parameter :: &
       ncoef_got=10,          &      ! # coefficient for properties
       nexap_got=10                  ! # exact solution parameters
  character(150)           :: &
       fil_rstar_got
  real(rp),      parameter :: &
       zensi = epsilon(1.0_rp)
!--BEGIN REA GROUP
  !------------------------------------------------------------------------
  ! Physical problem: read in got_reaphy
  !------------------------------------------------------------------------
  integer(ip)              :: &
       kfl_diffu_got,         &      ! Existence diffusion
       kfl_difun_got,         &      ! Diffusion function
       kfl_forme_got,         &      ! Conservative/non-conservative
       kfl_probl_got,         &      ! Problem to be solved
       kfl_timei_got,         &      ! Existence of du/dt of da/dt
       kfl_timec_got,         &      ! Existence of da/dt in continuity
       kfl_timem_got,         &      ! Existence of du/dt in momentum
       kfl_velfu_got                 ! Velocity function
  real(rp) ::&
       ddrop_got,             &      ! Droplet diameter (d)
       deair_got,             &      ! Air density (rho_a)
       densi_got,             &      ! Density density (rho)
       diffu_got(3),          &      ! Diffusion value (k)
       gravi_got(3),          &      ! Gravity vector
       grnor_got,             &      ! Gravity norm
       leinf_got,             &      ! Characteristic length
       muair_got,             &      ! Air viscosity (mu_a)
       veair_got                     ! Characteristic air velocity

  !------------------------------------------------------------------------
  ! Numerical problem: read in got_reanut
  !------------------------------------------------------------------------
  integer(ip)              :: &
       kfl_algor_got,         &      ! Algorithm (monolithic/Block Gauss-Seidel9
       kfl_artif_got,         &      ! Artificial viscosity
       kfl_coupl_got,         &      ! Equations coupling
       kfl_dttyp_got,         &      ! Local strategy of time step
       kfl_ellen_got,         &      ! =0,1 for min/max element length
       kfl_linea_got,         &      ! Linearization (RHS=0, Picard=1, Newton=2)
       kfl_normc_got,         &      ! Norm of convergence
       kfl_penal_got,         &      ! Momentum penalization
       kfl_sgsco_got,         &      ! Stabilization convection tracking
       kfl_sgsti_got,         &      ! Stabilization time tracking
       kfl_shocc_got,         &      ! Shock capturing continuity
       kfl_shocm_got,         &      ! Shock capturing momentum
       kfl_staty_got,         &      ! Stabilization type (SUPG/ASGS)
       kfl_taust_got,         &      ! Tau strategy
       kfl_tiacc_got,         &      ! Temporal accuracy
       kfl_weigc_got,         &      ! Weight continuity residual 
       kfl_weigm_got,         &      ! Weight momentum residual 
       itart_got,             &      ! Stop artificial viscosity
       itshc_got,             &      ! Stop shock capturing continuity
       itshm_got,             &      ! Stop shock capturing momentum
       mibgs_got,             &      ! Iteration block GS
       miinn_got,             &      ! Max # of iterations
       misgs_got,             &      ! Max # of SGS iterations
       neule_got,             &      ! # Euler iterations
       npica_got
  real(rp)                 :: &
       artif_got(2),          &      ! Artificial viscosity coefficient
       cotol_got,             &      ! Convergence tolerance
       cutof_got,             &      ! Alpha relative cut-off 
       penal_got,             &      ! Penalization factor
       relax_got,             &      ! Relaxation parameter
       relgs_got,             &      ! Relaxation BGS
       relsg_got,             &      ! Relaxation parameter of subgrid scale
       safet_got,             &      ! Safety factor for time step
       shock_got,             &      ! Shock capturing parameter
       sstol_got,             &      ! Steady state tolerance
       staco_got(3),          &      ! Stability constants
       tobgs_got,             &      ! Tolerance BGS
       tosgs_got,             &      ! Subgrid scale tolerance
       xmaxi_got(3),          &      ! Max bounding box
       xmini_got(3)                  ! Min bounding box

  !------------------------------------------------------------------------
  ! Boundary conditions: read in got_reabcs
  !------------------------------------------------------------------------
  integer(ip), pointer     :: &
       kfl_fixno_got(:)              ! Nodal fixity 
  real(rp),    pointer     :: &
       bvess_got(:,:)                ! Essential bc values

  !------------------------------------------------------------------------
  ! Output and Postprocess: read in got_reaous
  !------------------------------------------------------------------------
  integer(ip)              ::    &
       kfl_exacs_got,            &    ! Exact solution
       npp_bound_got                  ! Postprocess boundary conditions
  real(rp) ::&
       pos_cutof_got                  ! Cut off of postprocess
!--END REA GROUP
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  real(rp) ::                 &
       expar_got(nexap_got)          ! Exact solution parameters

  !
  ! Dimensions
  !
  integer(ip)              :: &
       ndofn_got(3),          &      ! # of d.o.f. of the NSI problem
       ndof2_got(3),          &      ! ndofn_got*ndofn_got
       ncomp_got,             &      ! Number of components of the droplet velocity (NSI)
       nunkn_got(2),          &      ! # of unknonws ndofn*npoin  
       nevat_got,             &      ! Element matrix dim.=(ndime+1)*mnode
       nzmat_got(2),          &      ! Matrix size
       nzrhs_got(2),          &      ! RHS size
       nzpre_got(2)                  ! Preconditioner size
  !
  ! Internal variables
  !
  integer(ip)              :: &
       ittot_got,             &      ! Total number of iteration
       kfl_goite_got,         &      ! Keep iterating
       kfl_resid_got,         &      ! If droplet velocity residual is required for post.
       kfl_tiaor_got,         &      ! Original time accuracy
       kfl_stead_got,         &      ! Steady-state has been reached 
       ivari_got                     ! Equation being solved (momentum and/or continuity) 
  integer(ip), allocatable :: &
       itsta_got(:)                  ! Statistics sgs
  real(rp)                 :: &
       dtcri_got,             &      ! Critical time step
       dtinv_got,             &      ! 1/dt,
       kfact_got,             &      ! K factor (K)
       resiv_got,             &      ! Velocity residual
       resic_got,             &      ! Water volume fraction
       resgs_got,             &      ! Subgrid scale residual
       exvdr_got(2),          &      ! Exact solution: u
       excdr_got,             &      ! Exact solution: alpha
       exgvd_got(2,2),        &      ! Exact solution: grad(u)
       exgcd_got(2),          &      ! Exact solution: grad(alpha)
       almax_got,             &      ! Maximum alpha  
       vemax_got,             &      ! Maximum Velocity  
       tamin_got,             &      ! Minimum tau
       tamax_got,             &      ! Maximum tau
       dimin_got,             &      ! Minimum diffusion
       dimax_got,             &      ! Maximum diffusion
       timin_got,             &      ! Minimum time step
       timax_got                     ! Maximum time step
  real(rp),    pointer     :: &
       veloc_got(:,:),        &      ! Velocity read from file
       cdold_got(:),          &      ! For residual
       vdold_got(:,:),        &      ! For residual
       diffm_got(:),          &      ! Diffusion for momentum equation
       elcod_got(:,:),        &      ! ELCOD: Elemental operations
       elvel_got(:,:),        &      ! ELVEL: Elemental operations
       elvdr_got(:,:,:),      &      ! ELVDR: Elemental operations
       elcdr_got(:,:),        &      ! ELCDR: Elemental operations
       eldif_got(:)                ! ELDIF: Elemental operations

  integer(ip), target      :: &
       itsol_got                     ! Solver iterations
  real(rp),    allocatable :: &    
       resis_got(:,:)                ! Subgrid scale inner residual
  !
  ! Subgrid scales of primary variables
  !
  type(r2p),   pointer     :: &
       vdsgs_got(:)                  ! Velocity subgrid scale
  type(r1p),   pointer     :: &
       cdsgs_got(:)                  ! Water volume fraction subgrid scale
  !
  ! Solver
  ! 
  real(rp),    pointer     :: &
       gslhs_got(:),          &      ! Gauss-Seidel: LHS
       gsrhs_got(:),          &      ! Gauss-Seidel: RHS
       gsunk_got(:),          &      ! Gauss-Seidel: Unknown
       gspre_got(:),          &      ! Gauss-Seidel: preconditioner
       unbgs_got(:)                  ! Gauss-Seidel: unknown
  real(rp),    pointer     :: &
       amatr_got(:)                  ! Linear matrix

end module def_gotita
