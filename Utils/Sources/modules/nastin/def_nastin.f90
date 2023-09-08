module def_nastin
  !------------------------------------------------------------------------
  !    
  ! Heading for the incompressible NASTIN routines
  !
  !------------------------------------------------------------------------
  use def_kintyp
  use def_coupli,        only :  typ_color_coupling

  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  integer(ip),   parameter :: &
       NSI_INCREMENTAL_PROJECTION          =  2, &
       NSI_PREDICTOR_CORRECTOR             =  3, &
       NSI_BLOCK_GAUSS_SEIDEL              =  4, &
       NSI_MOMENTUM                        =  2, &
       NSI_CONTINUITY                      =  3, &
       NSI_MOMENTUM_AND_CONTINUITY         =  1, &
       NSI_INCOMPRESSIBLE                  =  0, &
       NSI_COMPRESSIBLE                    =  1, &
       NSI_LOW_MACH                        =  3, &
       NSI_ANALYTICAL_HYDROSTATIC_PRESSURE =  1, &
       NSI_PDE_HYDROSTATIC_PRESSURE        =  2, &
       NSI_GALERKIN                        = -1, &
       NSI_ASGS                            =  0, &
       NSI_OSS                             =  1, &
       NSI_SPLIT_OSS                       =  2, &
       NSI_ALGEBRAIC_SPLIT_OSS             =  3, &
       NSI_CONVECTION_NON_CONSERVATIVE     =  0, &
       NSI_CONVECTION_CONSERVATIVE         =  1, &
       NSI_CONVECTION_SKEW                 =  2, &
       NSI_CONVECTION_EMAC                 =  3, &
       NSI_DIRICHLET_ELEMENT               =  0, &
       NSI_DIRICHLET_MATRIX                =  1, &
       NSI_DIRICHLET_ALGORITHM             =  2, &
       NSI_GLOBAL_TO_LOCAL                 =  1,&
       NSI_LOCAL_TO_GLOBAL                 =  2,&
       NSI_LUMPED_MASS                     =  0,&
       NSI_CONSISTENT_MASS                 =  1
  
  integer(ip),   parameter :: &
       lun_bound_nsi = 110, lun_stasg_nsi = 112, lun_cvgsg_nsi = 113, &
       lun_psmat_nsi = 114, lun_refer_nsi = 117, &
       lun_recvg_nsi = 118, lun_lmach_nsi = 120, lun_dynin_nsi = 122, &
       lun_dynou_nsi = 123, lun_dynlo_nsi = 124, lun_dynre_nsi = 125, &
       lun_windt_nsi = 1100
  character(150)           :: &
       fil_conve_nsi,         &      ! Convergence file name
       fil_dynin_nsi,         &      ! Input dynamic model
       fil_dynou_nsi,         &      ! Output dynamic model
       fil_windt_nsi                 ! Wind turbine output file 
  integer(ip),   parameter :: &
       nvars_nsi=25,          &      ! # set variables
       nvart_nsi=10,          &      ! # times for postprocess
       nvarw_nsi=10,          &      ! # witness point
       nvarp_nsi=40,          &      ! # postprocess variables
       ncoef_nsi=10,          &      ! # coefficient for properties
       mtabl_nsi=500,          &      ! maximum number of tabulated ct cp coefs for disk properties
       mforc_material_nsi=21, &      ! Maximum number of material force parameters 
       mflow_nsi=500,          &     ! # max number of flow rate conditions
       ivari_nsi_mom=1,       &      ! Momentum 
       ivari_nsi_cont=2,      &      ! Continuity 
       ivari_nsi_corre=4,     &      ! end of step correction 
       ivari_nsi_hydro=5,     &      ! hydrostatic state
       ivari_nsi_divfree=6,   &      ! Divergence Free
       mmsgs_nsi=50                  ! Max. number of iterations for SGS convergence statistics
  
  character(150)           :: &
       fil_rstar_nsi
  real(rp),      parameter :: &
       zensi = epsilon(1.0_rp)

  !--BEGIN REA GROUP
  !------------------------------------------------------------------------
  !
  ! Physical problem: read in nsi_reaphy
  !
  !------------------------------------------------------------------------

  integer(ip)               :: &
       kfl_advec_nsi,          &      ! Existence of (u.grad)u
       kfl_convection_type_nsi,&      ! Type of convection
       kfl_colev_nsi,          &      ! Coupling with LEVELS
       kfl_cotem_nsi,          &      ! Coupling with TEMPER
       kfl_cotur_nsi,          &      ! Coupling with TURBUL
       kfl_timei_nsi,          &      ! Existence of du/dt
       kfl_fvfua_nsi,          &      ! Frame angular velocity function  
       kfl_fvful_nsi,          &      ! Frame linear velocity function  
       kfl_grtur_nsi,          &      ! Add grad(K) to momemtum equations
       kfl_regim_nsi,          &      ! Flow regime: incompressible/compressible 
       kfl_dynco_nsi,          &      ! Dynamical coupling
       kfl_visco_nsi,          &      ! Viscous term 
       kfl_prthe_nsi,          &      ! Thermodynamic pressure calculation
       kfl_surte_nsi,          &      ! Include surface tension
       kfl_force_nsi,          &      ! Force term of the momentum equations
       kfl_mfrco_nsi,          &      ! Mass flow rate control activation flag
       kfl_bnods_nsi,          &      ! boundary nodes defined
       kfl_hydro_gravity_nsi,  &      ! Add hydrostatic gravity to NS
       kfl_fscon_nsi,          &      ! FS consistent algorithm (solve a Poisson each step)
       mfrse_nsi,              &      ! Set from which the mass flow rate is calculated
       nbval_nsi,              &      ! number of boundary values for time-space boundary from file
       nbtim_nsi,              &      ! number of time instances for time-space boundary from file
       nfiel_nsi(2),            &      ! Fields assignement 
       nbnod_nsi                      ! Number of nodes on boundary

  real(rp)                  :: &
       boube_nsi    ,          &      ! Boussinesq volume expansion
       bougr_nsi    ,          &      ! Boussinesq gravity
       boutr_nsi    ,          &      ! Boussinesq reference temperature
       lowtr_nsi    ,          &      ! Low Mach reference temperature
       facca_nsi(3) ,          &      ! Frame angular acceleration vector       
       faccl_nsi(3) ,          &      ! Frame linear acceleration vector  
       fadia_nsi(3) ,          &      ! Frame angular acceleration direction 
       fadil_nsi(3) ,          &      ! Frame linear acceleration direction 
       fvdil_nsi(3) ,          &      ! Frame linear velocity direction     
       fanoa_nsi    ,          &      ! Frame angular acceleration norm   
       fanol_nsi    ,          &      ! Frame linear acceleration norm      
       fcons_nsi    ,          &      ! Convection term
       frotc_nsi(3) ,          &      ! Frame rotation center
       centr_nsi    ,          &      ! Centrifugal force
       fvdia_nsi(3) ,          &      ! Frame angular velocity direction 
       fvela_nsi(3) ,          &      ! Frame angular velocity vector    
       fvell_nsi(3) ,          &      ! Frame linear velocity vector        
       fvins_nsi    ,          &      ! Viscous term
       fvnoa_nsi    ,          &      ! Frame angular velocity norm     
       fvnol_nsi    ,          &      ! Frame linear velocity norm 
       fvpaa_nsi(6) ,          &      ! Frame angular velocity parameters  
       fvpal_nsi(6) ,          &      ! Frame linear velocity parameters  
       gravi_nsi(3) ,          &      ! Gravity vector
       gravb_nsi(3) ,          &      ! Gravity vector for Boussinesq coupling
       grnor_nsi    ,          &      ! Gravity norm
       turbu_nsi(2) ,          &      ! Turbulence parameters
       heihy_nsi    ,          &      ! Height for hydrostatic pressure
       surte_nsi    ,          &      ! Surface tension coeficient (sigma)
       mfrub_nsi    ,          &      ! Target bulk velocity when mass flow control activated
       nbtdt_nsi,              &      ! Time-step of the turbulent inlet database
       fsrot_nsi,              &      ! FS rotational factor
       ubpre_nsi    ,          &      ! Bulk velocity from previous time-step
       mfccf_nsi                      ! Coefficient for the mass flow control formula

  integer(ip), pointer      :: &
       lforc_material_nsi(:)  ,&      ! List of material force 
       ntabl_nsi(:)           ,&      ! number of tabulated ct and cp parameters (wind turbines)
       ntabr_nsi(:)                   ! number of tabulated rotational parameters (wind turbines)
  real(rp),    pointer      :: &
       xforc_material_nsi(:,:),&        ! Material force parameters
       velta_nsi(:,:)         ,&        ! tabulated input velocity 
       thrta_nsi(:,:)         ,&        ! tabulated thrust coeff
       powta_nsi(:,:)         ,&        ! tabulated power  coeff
       veave_nsi(:,:)         ,&        ! averaged velocity
       radiu_nsi(:,:)         ,&        ! tabulated dimensional radius
       forcn_nsi(:,:)         ,&        ! tabulated normal force distribution
       forct_nsi(:,:)                   ! tabulated tangential force distribution

  real(rp),    pointer         :: &
       bntab_nsi(:,:),            &     ! boundary nodes table for time-space boundary from file
       bnval_nsi(:,:)                   ! boundary values table for time-space boundary from file

  integer(ip),    pointer      :: &
       iboun_nsi(:)                     ! boundary correspondence for time-space boundary from file


  !
  ! Fluid properties
  !
  real(rp)                  :: &
       sphea_nsi,              &      ! Specific heat (Cp)
       prthe_nsi,              &      ! Thermodynamics pressure (cst or initial)
       tmass_nsi                      ! Initial mean density related to initial mass(low Mach)

  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in nsi_reanut
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_penal_nsi,         &      ! Penalization
       kfl_prepe_nsi,         &      ! Pressure Penalization - Used to avoid temporal pressure oscilations with local dt
       kfl_dttyp_nsi,         &      ! Local strategy of time step
       kfl_ellen_nsi,         &      ! =0,1 for min/max element length
       kfl_relax_nsi,         &      ! Velocity relaxation strategy
       kfl_relap_nsi,         &      ! Pressure relaxation strategy
       kfl_sgsco_nsi,         &      ! Stabilization convection tracking
       kfl_sgsti_nsi,         &      ! Stabilization time tracking
       kfl_sgsac_nsi,         &      ! Stabilization time tracking accuracy
       kfl_sgsli_nsi,         &      ! Stabilization tracking convection linearization
       kfl_sgscp_nsi,         &      ! Coupling of SGS with grid scale
       kfl_shock_nsi,         &      ! Shock capturing type 
       kfl_tiacc_nsi,         &      ! Temporal accuracy
       kfl_normc_nsi,         &      ! Norm of convergence
       kfl_refer_nsi,         &      ! Difference between solutions
       kfl_linea_nsi,         &      ! Linearization (RHS=0, Picard=1, Newton=2, ExactNewton(adj)=3)
       kfl_tisch_nsi,         &      ! Time integration scheme
       kfl_algor_nsi,         &      ! Type of algorithm ( Monolithic, gauss-seidel, etc.)
       kfl_predi_nsi,         &      ! Predictor corrector
       kfl_taush_nsi,         &      ! Schur complement. Tau strategy 
       kfl_ellsh_nsi,         &      ! Schur complement. =0,1 for min/max element length 
       kfl_updpr_nsi,         &      ! Pressure update
       kfl_intpr_nsi,         &      ! Treatment of the pressure term
       kfl_assem_nsi,         &      ! Assembly type (1,2)
       kfl_asbou_nsi,         &      ! Assembly type (1,2)
       kfl_taust_nsi,         &      ! Tau strategy
       kfl_stabi_nsi,         &      ! Orthogonal SGS
       kfl_limit_nsi,         &      ! Limiter for split OSS
       kfl_trres_nsi,         &      ! Transient residual
       kfl_prtre_nsi,         &      ! Pressure treatment (explicit/implicit)
       kfl_matdi_nsi,         &      ! Dirichlet bc on matrix
       kfl_intfo_nsi,         &      ! Internal force calculation (=0: integral, 1=residual based)
       kfl_press_nsi,         &      ! Integrate pressure term in momentum equations
       kfl_bubbl_nsi,         &      ! Pressure bubble
       kfl_stain_nsi,         &      ! Time step to start inner iterations
       kfl_immer_nsi,         &      ! Immersed boundary method
       kfl_grad_div_nsi,      &      ! Use grad an div matrices, also Laplacian despite the name does not reflect it
       misgs_nsi,             &      ! Max # of SGS iterations
       npica_nsi,             &      ! Number of Picard iteration (Newton's lin.)
       neule_nsi,             &      ! # of Euler time steps
       kfl_meshi_nsi,         &      ! Mesh interpolator activation flag
       kfl_savco_nsi,         &      ! Save linear matrix
       kfl_corre_nsi,         &      ! Fractional step correction-like 
       kfl_sosch_nsi,         &      ! Schur complement solver
       kfl_modfi_nsi,         &      ! Modify kfl_fixno_nsi
       kfl_expco_nsi,         &      ! Treat the convective term explicitly, that is, assemble the matrix only in the first ortomin iteration
       kfl_addpr_nsi,         &      ! Add contribution due to pressure in matrix side (do nothing BC') on wall law boundaries        
       kfl_grvir_nsi,         &      ! Add  viscous gradient contribution inside the residual
       kfl_wlare_nsi,         &      ! Initialize time-averaged velocity after restart
       kfl_hydro_nsi,         &      ! Hydrostatic initial state
       kfl_update_hydro_nsi,  &      ! When to update hydrostatic pressure
       kfl_hydro_interface_nsi, &    ! Interface height computation
       mitri_nsi,             &      ! Maximum number of Richardson iterations
       kfl_adres_nsi,         &      ! FS, Tau solver: adaptive residual 
       kfl_incre_nsi,         &      ! Solve pressure in incremental form
       kfl_nota1_nsi,         &      ! Does not stabilize of convective reactive and coriolis term in momentum eq.
       kfl_ini_ts_guess_order_nsi,&  ! Order of the initial guess extrapolation at each time step  - for the moment only velocity. 
       kfl_vector_nsi,        &      ! Vectorized version of Alya
       kfl_press_stab_nsi,    &      ! Pressure stabilization: dt(0) or tau(1)
       kfl_stop_by_wit_nsi,   &      ! Stop by convergence of witness points - uses velocity
       kfl_massm_nsi                 ! Consistent mass matrix
  
  real(rp)                 :: &
       penal_nsi,             &      ! Penalization factor
       prepe_nsi,             &      ! Pressure penalization factor
       dtcri_nsi,             &      ! Critical time step
       staco_nsi(4),          &      ! Stability constants
       shock_nsi,             &      ! Shock capturing parameter
       safet_nsi,             &      ! Safety factor for time step
       bemol_nsi,             &      ! Integration of convective term by parts
       sstol_nsi,             &      ! Steady state tolerance
       cotol_nsi,             &      ! Convergence tolerance
       resid_nsi,             &      ! Residual for outer iterations (u)
       resip_nsi,             &      ! Residual for outer iterations (p)
       weigh_nsi,             &      ! Weight of dU/dt in the residual
       relax_nsi,             &      ! Relaxation parameter velocity
       relap_nsi,             &      ! Relaxation parameter pressure
       relsg_nsi,             &      ! Relaxation parameter of subgrid scale
       tosgs_nsi,             &      ! Subgrid scale tolerance
       strec_nsi,             &      ! Adaptive dt: Stretching factor
       dampi_nsi,             &      ! Adaptive dt: damping
       epsht_nsi,             &      ! Adaptive dt: eps_R
       epstr_nsi,             &      ! Adaptive dt: eps_A
       xfree_nsi,             &      ! X Coordinate of the plane where to free
       safex_nsi,             &      ! Time function parameter for safety factor
       adres_nsi,             &      ! FS, Tau solver: adaptive residual factor
       toler_nsi,             &      ! FS, Tau solver: adaptive residual factor
       safma_nsi,             &      ! Maximum safety factor
       safeo_nsi,             &      ! Initial safety factor
       saflo_nsi,             &      ! Minimum global safety factor for local time steps
       gamma_nsi                     ! Gamma factor for pressure in momentum equation
 
  !------------------------------------------------------------------------
  !
  ! Boundary conditions: read in nsi_reabcs
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_confi_nsi,         &       ! Confined flow
       kfl_local_nsi,         &       ! Local system of reference
       kfl_conbc_nsi,         &       ! Constant b.c.
       kfl_syntu_nsi,         &       ! Synthetic eddy method (SEM)
       kfl_initi_nsi,         &       ! Initial 
       kfl_inico_nsi,         &       ! Solve coarse grid system
       kfl_inipr_nsi,         &       ! Initial pressure
       kfl_nopen_nsi,         &       ! No penetration condition
       kfl_cadan_nsi,         &       ! Coupling with ADAN
       kfl_aiobo_nsi,         &       ! Alya IO boundary to couple with ADAN 
       itebc_nsi,             &       ! Initial step for boundary condition
       nodpr_nsi,             &       ! Node on which pressure is prescribed
       exfpr_nsi,             &       ! Extends fixpr one layer of elements
       neddy_nsi,             &       ! Number of eddies in the inlet box for SEM
       kfl_imppr_nsi,         &       ! Imposses pressure in nodes w/ fixpr>0 
       kfl_divcorrec_nsi              ! Correction div(u)=0 to initial solution       
  real(rp) ::&
       delta_nsi,             &      ! Distance to the wall
       relbc_nsi,             &      ! Boundary condition relaxation
       valpr_nsi,             &      ! Pressure value 
       velin_nsi(3),          &      ! Initial constant velocity
       poise_nsi(6),          &      ! Parameters for Poiseuille law
       divcorrec_nsi                 ! div(u)=0 parameter (alpha)
  real(rp) ::&
       press_cadan_nsi               ! Pressure coming from ADAN (ADAN COUPLING)
  real(rp), allocatable    :: &
       Q_cadan_nsi(:),        &      ! Flow to send to ADAN
       P_cadan_nsi(:)                ! Pressure to send to ADAN
  type(bc_nodes), pointer  :: &     
       tncod_nsi(:)                  ! Node code type
  type(bc_nodes), pointer  :: &     
       tgcod_nsi(:)                  ! Geometrical node code type
  type(bc_bound), pointer  :: &     
       tbcod_nsi(:)                  ! Boundary code type

  integer(ip), pointer  ::&
       kfl_flow_rate_codes_nsi(:),& ! Flow rates codes 
       kfl_flow_rate_normal_nsi(:),& ! Flow rates normals 
       kfl_flow_rate_stfn_nsi(:)    ! Flow rates space time function number       
  real(rp), pointer  ::&
       flow_rate_values_nsi(:),&     ! Flow rate values
       flow_rate_normals_nsi(:,:)     ! Flow rate values

  !------------------------------------------------------------------------
  !
  ! Optimization problem
  !
  !------------------------------------------------------------------------

  real(rp),   pointer ::                   &
       resdiff_nsi(:,:),                   & ! Partial Derivative of residual w.r.t design var
       dcost_dx_nsi(:),                    & ! Partial Derivative of objective function w.r.t coordinates
       costdiff_nsi(:)                       ! Partial Derivative of objective function w.r.t design var
  
  !------------------------------------------------------------------------
  !
  ! Output and Postprocess: read in nsi_reaous
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_exacs_nsi,         &      ! Exact solution for the NS eq.
       kfl_exfix_nsi,         &      ! Fixity imposed for exact solution
       kfl_psmat_nsi,         &      ! Postscript of matrix profile
       kfl_inert_nsi                 ! Velocity in inertial frame of ref.
  real(rp) ::&
       expar_nsi(10),           &    ! Exact solution parameters
       cloth_nsi,               &    ! CLO
       metab_nsi,               &    ! MET
       wetme_nsi,               &    ! WME
       ambie_nsi,               &    ! TA
       radia_nsi,               &    ! TR
       relat_nsi,               &    ! RH
       avtim_nsi,               &    ! Start averaging time
       entim_nsi                     ! End ensemble time
  !--END REA GROUP
  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------
  !
  ! Boundary conditions
  !
  integer(ip), pointer     :: &
       kfl_fixno_nsi(:,:),    &      ! Nodal fixity 
       kfl_fixpr_nsi(:,:),    &      ! Nodal fixity for the pressure Schur complement
       kfl_fixpp_nsi(:,:),    &      ! Nodal fixity for the pressure 
       kfl_fixbo_nsi(:),      &      ! Element boundary fixity
       kfl_fixrs_nsi(:),      &      ! Reference system for the BV
       kfl_funno_nsi(:),      &      ! Functions for node bc
       kfl_funbo_nsi(:),      &      ! Functions for node bc
       kfl_wlawf_nsi(:),      &      ! Flag to identify if a point is part of a wall law boundary
       lexlo_nsi(:,:),        &      ! List for the exchange location for wall law
       ielem_wel(:),          &
       lbpse(:)                      ! List of boundary sets passed to nodes
  real(rp),    pointer     :: &
       bvess_nsi(:,:,:),      &      ! Essential velocity bc values
       bpess_nsi(:,:,:),      &      ! Essential bc values for pressure
       bvnat_nsi(:,:,:),      &      ! Natural bc values
       skcos_nsi(:,:,:),      &      ! Cosine matrices of NS skew systems
       velel_nsi(:,:),        &      ! Velocity at the exchange location for wall law  
       massb_nsi(:),          &      ! Boundary surface
       notra_nsi(:,:),        &      ! Traction on the boundary nodes for postprocess
       avntr_nsi(:,:),        &      ! Time average traction on the boundary nodes - variational calculation 
       avgtr_nsi(:,:),        &      ! Time average traction on the boundary nodes - using velocity gradients
       shape_wel(:,:),        &      ! shape functions associated to excange location (for implicit)
       btrac_nsi(:,:),        &      ! Traction on the boundary nodes from auxiliary RANS simulation
       tracr_nsi(:,:),        &      ! Traction calculated in auxiliary RANS simulation
       tluav_nsi(:,:)                ! Average velocity for two-layer coupling

  type(typ_color_coupling) :: wallcoupling
  !
  ! Dimensions
  !
  integer(ip)              :: &
       ndofn_nsi(2),          &      ! # of d.o.f. of the NSI problem
       ndof2_nsi(2),          &      ! ndofn_nsi*ndofn_nsi
       ncomp_nsi,             &      ! Number of components of the velocity (NSI)
       nprev_nsi,             &      ! Previous time step or iteration
       nunkn_nsi(2),          &      ! # of unknonws ndofn*npoin  
       nevat_nsi,             &      ! Element matrix dim.=(ndime+1)*mnode
       nzsol_nsi,             &      ! Matrix size (per d.o.f.)
       nzmat_nsi(2),          &      ! Matrix size
       nzrhs_nsi(2),          &      ! RHS size
       nzpre_nsi(2),          &      ! Preconditioner size
       lperp_nsi(8),          &      ! List of periodic prescribed pressure 
       kfl_perip_nsi,         &      ! If pressure is prescribed on periodic nodes
       kfl_dodem_nsi
  !
  ! Internal variables
  !
  logical(lg)              :: &
       NSI_MONOLITHIC,        &      ! Monolithic algorithm
       NSI_SCHUR_COMPLEMENT,  &      ! Schur complement algorithm
       NSI_FRACTIONAL_STEP
  integer(ip)              :: &
       ittot_nsi,             &      ! Total number of iteration
       kfl_resid_nsi,         &      ! If velocity residual is required for post.
       kfl_grvis_nsi,         &      ! If velocity gradients exist
       kfl_goite_nsi,         &      ! Keep iterating
       kfl_rmom2_nsi,         &      ! Off-diagonal part of momentum operator exists
       kfl_p1ve2_nsi,         &      ! Off-diagonal part of momentum test function exists
       ndbgs_nsi,             &      ! Number dof for BGS
       kfl_stead_nsi,         &      ! Steady-state has been reached 
       kfl_tiaor_nsi,         &      ! Original time accuracy
       kfl_sgste_nsi,         &      ! Temperature subgrid scale considered
       kfl_autom_nsi,         &      ! Automatic boundaries
       ivari_nsi,             &      ! Equation being solved (momentum and/or continuity)
       iteqn_nsi(2),          &      ! Internal iterations for momentum+continuity
       itbgs_nsi,             &      ! BGS iteration number
       nbdfp_nsi,             &      ! Number of terms in the temporal derivative
       kfl_exist_fixi7_nsi,   &      ! exists nodal fixity of type 7 
       kfl_exist_fib20_nsi,   &      ! exists boundary fixity of type 20 
       itsta_nsi(mmsgs_nsi)          ! Statistics sgs
  real(rp)                 :: &
       dtinv_nsi,             &      ! 1/dt , from vers 772 theta is now longer included in dtinv_nsi
       dtsgs_nsi,             &      ! 1/(theta'*dt)
       err01_nsi(2),          &      ! L1 error u
       err02_nsi(2),          &      ! L2 error u
       err0i_nsi(2),          &      ! Linf error u
       err11_nsi(2),          &      ! L1 error grad(u)
       err12_nsi(2),          &      ! L2 error grad(u)
       err1i_nsi(2),          &      ! Linf error grad(u)
       corio_nsi,             &      ! Coriolis force
       pabdf_nsi(10),         &      ! BDF parameters, actually now (vers 772) we will extend it also for CN 
       rgsve_nsi,             &      ! residual BGS velocity
       rgspr_nsi,             &      ! residual BGS pressure
       resin_nsi(2),          &      ! Algebraic inner residual
       resou_nsi(2),          &      ! Algebraic outer residual
       resss_nsi(2),          &      ! Algebraic steady state residual
       reinf_nsi(2),          &      ! Algebraic Linf residual 
       tamin_nsi,             &      ! Min tau
       tamax_nsi,             &      ! Max tau
       vemin_nsi,             &      ! Min velocity
       vemax_nsi,             &      ! Max velocity
       prmin_nsi,             &      ! Min pressure
       prmax_nsi,             &      ! Max pressure
       pcoef_nsi,             &      ! Pressure coefficient 1-R/Cp
       relpa_nsi(2),          &      ! Relaxation parameter
       cputi_nsi(10),         &      ! CPU time
       cpu_ass_sol_nsi(4),    &      ! CPU time assembly and solver at each iteration
       cputi_assembly_nsi(10),&      ! COU time for element assembly
       gamth_nsi,             &      ! gamma=Cp/(Cp-R)
       xmass_nsi,             &      ! Low-Mach: Mass computed from state equation
       actav_nsi,             &      ! Accumulated time for averaging
       difve_nsi,             &      ! Velocity residual w/r reference solution
       difpr_nsi,             &      ! Pressure residual w/r reference solution
       vinvt_nsi(4),          &      ! for lowmac, = integ(1/T)
       hydro_nsi,             &      ! Hydrostatic z-plane
       dtmax_nsi,             &      ! for local time step stores the maximum time step
       rmsgs_nsi,             &      ! Maximum subgrid scale residual
       resgs_nsi(2),          &      ! Subgrid scale residual (numerator/denominator)
       resis_nsi(2,mmsgs_nsi)        ! Subgrid scale inner residual
  real(rp),    pointer     :: &    
       veold_nsi(:,:),        &      ! Velocity for residual
       gradv_nsi(:,:),        &      ! velocity gradient (postprocess)
       unk2n_nsi(:,:),        &      ! Nastin second variables (pressure or density)
       dunkn_nsi(:),          &      ! Delta velocity
       dunkp_nsi(:),          &      ! Delta pressure
       avvel_nsi(:,:),        &      ! Average velocity
       avve2_nsi(:,:),        &      ! Average velocity**2
       avvxy_nsi(:,:),        &      ! Average vx*vy
       avpre_nsi(:),          &      ! Average pressure
       avpr2_nsi(:),          &      ! Average pressure**2
       avtan_nsi(:,:),        &      ! Average tangential force
       avmut_nsi(:),          &      ! Average turbulent viscosity
       avstx_nsi(:,:),        &      ! Average stress mu_t * grad(u)
       avsty_nsi(:,:),        &      ! Average stress mu_t * grad(v)
       avstz_nsi(:,:),        &      ! Average stress mu_t * grad(w)
       envel_nsi(:,:),        &      ! Ensemble velocity
       enve2_nsi(:,:),        &      ! Ensemble velocity**2
       envxy_nsi(:,:),        &      ! Ensemble vx*vy
       enpre_nsi(:),          &      ! Ensemble pressure
       enpr2_nsi(:),          &      ! Ensemble pressure**2
       entan_nsi(:,:),        &      ! Ensemble tangential force
       enmut_nsi(:),          &      ! Ensemble turbulent viscosity
       enstx_nsi(:,:),        &      ! Ensemble stress mu_t * grad(u)
       ensty_nsi(:,:),        &      ! Ensemble stress mu_t * grad(v)
       enstz_nsi(:,:),        &      ! Ensemble stress mu_t * grad(w)
       resch_nsi(:),          &      ! Schur complement residual
       remom_nsi(:,:),        &      ! Momentum residual
       prope_nsi(:,:),        &      ! Smoothed fluid Properties
       norle_nsi(:,:),        &      ! Normal to the zero Level Set 
       curle_nsi(:),          &      ! Curvature
       outflow_mass(:),       &      ! Outflow mass
       dt_rho_nsi(:),         &      ! Projection of dt / rho
       mu_rho_nsi(:),         &      ! Projection of dt / rho
       tau_nsi(:),            &      ! Projection of tau
       bubble_nsi(:),         &      ! Pressure bubble
       bubble_aqq_nsi(:),     &      ! Bubble matrix Aqq
       bubble_aqu_nsi(:,:),   &      ! Bubble matrix Aqu
       bubble_aqp_nsi(:,:),   &      ! Bubble matrix Aqp 
       bubble_bq_nsi(:),      &      ! Bubble RHS
       lagra_nsi(:,:,:),      &      ! Lagrange multiplier velocity
       tauib_nsi(:,:,:),      &      ! Lagrange multiplier tau
       vafor_nsi(:,:),        &      ! variational force
       avvaf_nsi(:,:),        &      ! time averaged variational force
       bupor_nsi(:,:)                ! forces due to porous media

  real(rp)                  ::       porfo_nsi(3)
  type(r3p),   pointer     :: &
       turmu_nsi(:)                  ! LES turbulent viscosity

  type nsimat
     integer(ip)           :: kfl_exist
     real(rp), pointer     :: Auu(:,:,:)
     real(rp), pointer     :: Aup(:,:)
     real(rp), pointer     :: bu(:)
  end type nsimat
  type(nsimat),   pointer  :: &
       intfo_nsi(:)                  ! Internal force
  !
  ! Solver
  ! 
  integer(ip), pointer     :: &
       kfl_fixno_div_nsi(:,:)        ! Nodal fixity for divergence free correction 
  integer(ip)              :: &
       nmauu_nsi,             &      ! Size of Auu
       nmaup_nsi,             &      ! Size of Aup
       nmapu_nsi,             &      ! Size of Apu
       nmapp_nsi,             &      ! Size of App
       poauu_nsi,             &      ! Pointer to Auu
       poaup_nsi,             &      ! Pointer to Aup
       poapu_nsi,             &      ! Pointer to Apu
       poapp_nsi,             &      ! Pointer to App
       nschu_nsi,             &      ! # Schur complement solves
       nmome_nsi                     ! # Momentum solves
  real(rp),    pointer     :: &
       Auu_nsi(:),            &      ! Auu
       Aup_nsi(:),            &      ! Aup
       Apu_nsi(:),            &      ! Apu
       App_nsi(:),            &      ! App
       amatr_nsi(:),          &      ! Linear matrix
       lapla_nsi(:),          &      ! Laplacian matrix
       cmama_nsi(:),          &      ! Consistent_mass_ matrix
       deltp_nsi(:),          &      ! Schur complement: Dp (used for mass correction)
       vepro_nsi(:,:),        &      ! Velocity projection
       grpro_nsi(:,:),        &      ! Pressure gradient projection
       prpro_nsi(:),          &      ! Pressure projection
       vepr2_nsi(:,:),        &      ! Velocity projection
       grpr2_nsi(:,:),        &      ! Pressure gradient projection
       prpr2_nsi(:)                  ! Pressure projection
  type(r1p),   pointer     :: &
       hydro_density_nsi(:)          ! Hydrostatic density

  integer(ip), pointer :: prout_nsi(:)
  
end module def_nastin
