module def_chemic
  !------------------------------------------------------------------------
  !****f* Partis/def_chemic
  ! NAME 
  !    def_chemic
  ! DESCRIPTION
  !    Heading for the Partis routines
  ! USES
  ! USED BY
  !    Almost all
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use mod_ADR,   only : ADR_typ

  implicit none
  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip), parameter ::                &
       lun_sized_chm = 1910,               &  ! Size distribution    
       lun_times_chm = 1911,               &  ! Time step info
       lun_time2_chm = 1912,               &  ! Time step targets
       lun_resu1_chm = 1921,               &  ! Result 1
       lun_resu2_chm = 1922,               &  ! Result 2
       lun_resu3_chm = 1923,               &  ! Result 3
       lun_remet_chm = 1930,               &  ! Meteo properties file
       lun_resou_chm = 1931,               &  ! Meteo source file
       lun_spcvg_chm = 1932                   ! Species convergence file
       
  character(150)                        :: &
       fil_sized_chm,                      &  ! Size distribution
       fil_times_chm,                      &  ! Time step    
       fil_time2_chm,                      &  ! Time step target
       fil_remet_chm,                      &  ! Meteo
       fil_resou_chm                          ! Meteo
  real(rp),      parameter :: &
       zepts = epsilon(1.0_rp)
  integer(ip),   parameter              :: &
       npara_chm=10,                       &   ! # parameters for sets
       npart_chm=10                            ! # temperature parameters
  logical(lg)                           :: &
       METEO_MODEL,                        &   ! Model
       DEFECT_EVOLUTION_MODEL                  ! Model
!--BEGIN REA GROUP
  !------------------------------------------------------------------------
  ! Physical problem: read in chm_reaphy
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_model_chm,                      & ! PDE-ODE model
       kfl_stagg_chm,                      & ! Staggered step combustion model
       kfl_timei_chm,                      & ! Existence of du/dt
       kfl_advec_chm,                      & ! Existence of (a.grad)u
       kfl_diffu_chm,                      & ! Existence of -div[k. grad(u)]
       kfl_react_chm,                      & ! Assemble reaction term on left hand side
       kfl_lhsas_chm,                      & ! Left hand side assembly of source terms
       kfl_activ_chm,                      & ! Activity coefficients present
       kfl_corve_chm,                      & ! Correction velocity
       kfl_norma_chm,                      & ! Normalize mass fractions
       kfl_sourc_chm,                      & ! Existence and type of source term
       kfl_meteo_chm,                      & ! Existence of meteo file
       kfl_ansou_chm,                      & ! Flag to identify the use of analytic sounding in Kessler
       kfl_tfles_chm,                      & ! Flag to activate the Thickened Flame model TFLES
       kfl_cotur_chm,                      & ! Flag to identify the RANS model for the CFI combustion model
       kfl_uncfi_chm,                      & ! Flag to activate unscale model for non-premixed flames in CFI combustion model
       kfl_tucfi_chm,                      & ! Flag to identify if variance is used for the CFI combustion model
       kfl_wallc_chm,                      & ! Flag to impose a zero source term at walls for the CFI combustion model
       kfl_radia_chm,                      & ! Flag to activate radiation model for CFI model
       kfl_field_chm(2),                   & ! Flag to activate the initialization by fields (1) Species fields, (2) Temperature field
       lawte_chm,                          & ! Law for temperature
       lawde_chm,                          & ! Law for density
       lawvt_chm,                          & ! Law for terminal velocity
       nclas_chm,                          & ! Number of particle classes
       nspec_chm,                          & ! Total number of species
       nodes_chm                             ! Number of ODE's
       
  real(rp)                              :: &
       sourc_chm,                          & ! Parameter for the source term
       radwt_chm,                          & ! Wall temprature for radiation model
       sorad_chm,                          & ! Source radius
       socen_chm(3),                       & ! Source center
       tfles_chm,                          & ! Thickening factor for TFLES
       flbet_chm                             ! beta for TFLES sensor

  integer(ip)                           :: &
       nreac_chm,                          & ! Number of reactions
       ncoef_chm,                          & ! Number of reaction coefficients
       kfl_arreh_chm,                      & ! Correction of Arrehnius coefficient with equivalence ratio
       stofu_chm(3)                          ! Which is fuel and which is oxygen for equivalence ratio

  real(rp)                              :: &
       denma_chm,                          & ! Material density
       radbi_chm,                          & ! Bi-molecular radius
       temma_chm(npart_chm),               & ! Material temperature
       boltz_chm,                          & ! Boltzmann constant
       strat_chm                             ! Stoichiometric mass ratio

  ! 
  integer(ip)                           :: &
       sponge_chm                            ! Viscosity amplification in sponge layer
  real(rp)                              :: &
     visco_factor_chm,                     & ! Amplification factor
     visco_axis(3),                        & ! Axis for amplification
     visco_range(2)                          ! Sponge range for amplification

  !
  ! Variables transmitted to nodes later
  integer(ip), pointer                  :: &
       lawdi_chm(:,:)                        ! Law diffusion
  real(rp),  pointer                    :: &      
       diffu_chm(:,:),                     & ! Diffusion constants
       radiu_chm(:),                       & ! Capture radius
       react_chm(:,:),                     & ! Reaction coefficients
       order_chm(:,:,:),                   & ! Reaction order (broken) parameters
       equil_chm(:,:),                     & ! 
       effic_chm(:,:),                     & ! Efficiency factor for reactions with chaperone
       stoic_chm(:,:,:),                   & ! Stoichiometric coefficients
       interaction_chm(:,:)                  ! Binary interaction coefficients


  type(i1p), pointer                    :: &
       lreac_chm(:)                          ! List of reactions

  real(rp), pointer                     :: &
       diame_chm(:),                       & ! METEO: Particle diameter
       rhopa_chm(:),                       & ! METEO: Particle Density
       shape_chm(:),                       & ! METEO: Particle Shape factor (Psi)
       spher_chm(:),                       & ! METEO: Particle Sphericity factor
       fract_chm(:)                          ! METEO: Particle Fraction
  !
  ! Note:
  ! Following variables are not send to slave nodes
  real(rp)                           :: & 
       afact_chm,                          & ! METEO:
       pbaro_chm,                          & ! METEO:
       cpcoe_chm,                          & ! METEO: 
       cvcoe_chm,                          & ! METEO: 
       cpliq_chm,                          & ! METEO:   
       adgam_chm,                          & ! METEO:`
       cpmli_chm,                          & ! METEO:   
       cpvap_chm,                          & ! METEO: 
       cvvap_chm,                          & ! METEO: 
       lhref_chm,                          & ! METEO: 
       teref_chm,                          & ! METEO: 
       rgasc_chm,                          & ! METEO: 
       rgava_chm                             ! METEO: 
  integer(ip)                           :: &
       kfl_benme_chm                        ! Meteo benchmark
  !
  real(rp),  pointer                    :: & 
       entha_chm(:,:),                     & ! Enthalpy for each species
       equiv_chm(:),                       & ! Equivalence ratio (mass)
       rspec_chm(:,:),                     & ! CFI model: species mass fraction for radiation model
       yscale_chm(:),                      & ! Scaled reaction progress to access database in partially premixed conditions
         cvar_chm(:),                      & ! Variance of the scaled reaction progress variable  
       flsen_chm(:),                       & ! Flame front sensor for THICKENED FLAME MODEL (TFLES)
       flspe_chm(:),                       & ! Flame speed for THICKENED FLAME MODEL (TFLES)
       flsgs_chm(:),                       & ! Subgrid scale wrinkling factor E for THICKENED FLAME MODEL (TFLES)
       flfac_chm(:),                       & ! Dyanmic thickened flame factor (DTFLES)
       flthi_chm(:)                          ! Flame thickness sensor for THICKENED FLAME MODEL (TFLES)

  character(5)                          :: &
       wprob_chm                             ! Problem name
  integer,parameter                     :: &
       maxsp_chm = 1000                       ! Max number of species tracked in combustion code


  !------------------------------------------------------------------------
  ! Numerical problem: read in chm_reanut
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_dttyp_chm,                      & ! Local time step strategy
       kfl_ellen_chm,                      & ! =0,1 for min/max element length
       kfl_taust_chm,                      & ! Tau calculation option
       kfl_shock_chm,                      & ! Shock capturing type 
       kfl_stabi_chm,                      & ! Stabilization strategy
       kfl_limit_chm,                      & ! Limiter
       kfl_tiacc_chm,                      & ! Temporal accuracy
       kfl_assem_chm,                      & ! Assembly strategy
       kfl_tibub_chm,                      & ! Time integration of bubble
       miinn_chm,                          & ! Max inner iterations
       neule_chm,                          & ! Number of Euler time steps 
       kfl_tisch_chm,                      & ! Time integration scheme
       kfl_normc_chm,                      & ! Norm of convergence
       kfl_coupl_chm,                      & ! Coupling of eqns
       kfl_dtcri_chm,                      & ! dt criteria
       kfl_dttar_chm,                      & ! dt target
       kfl_sgsti_chm,                      & ! Subscale time tracking
       kfl_negat_chm,                      & ! Startegy for negative concentrations
       kfl_posit_chm,                      & ! Startegy for too positive concentrations
       kfl_warni_chm,                      & ! Warn about points with zero sum mass
       kfl_meshi_chm,                      & ! Mesh interpolator activation flag
       kfl_temli_chm,                      & ! Flag to activate a T limiter to compute reaction rates
       kfl_gauss_chm,                      & ! Level of Gauss-Seidel update in combustion species
       kfl_spite_chm,                      & ! Intra-species iterations to improve shock capturing
       initial_fraction_step_chm,          & ! Fraction of total dt to evolve in staggered step scheme
       max_fixed_point_iterations_chm        ! Maximum iterations in fixed point for staggered step

  real(rp)                              :: &
       staco_chm(3),                       & ! Stability constants
       shock_chm,                          & ! Shock capturing parameter
       bemol_chm,                          & ! Bemol
       temli_chm,                          & ! Temperature limiter to compute reaction rates
       cotol_chm,                          & ! Convergence tolerance
       safet_chm,                          & ! Safety factor for time step
       chemical_time_factor,               & ! Safety factor exclusively for the source term
       cutof_chm,                          & ! Concentration cutoff for critical time computation
       sstol_chm,                          & ! Steady state tolerance
       strec_chm,                          & ! Adaptive dt: Stretching factor
       dampi_chm,                          & ! Adaptive dt: damping
       epsht_chm,                          & ! Adaptive dt: eps_R
       epstr_chm,                          & ! Adaptive dt: eps_A
       dtmin_chm,                          & ! Minimum time step
       dtmax_chm,                          & ! Maximum time step
       relax_chm,                          & ! Relaxation of update
       fixed_point_tolerance_chm,          & ! Tolerance for fixed point iterations in staggered step
       odeint_tolerance_chm,               & ! Tolerance for adaptive time stepping in ODE integration
       timestep_min_chm                      ! Minimum time step allowed in adaptive time stepping

  !------------------------------------------------------------------------
  ! Output and Postprocess: read in chm_reaous
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_sized_chm                         ! Size distribution
  integer(ip)                           :: &
       ipara_chm(npara_chm)                  ! Int parameters for sets
  real(rp)                              :: &
       rpara_chm(npara_chm),               & ! Real parameters for sets
       avtim_chm                             ! Accumulated time for time-averaging

  !------------------------------------------------------------------------
  ! Boundary conditions: read in chm_reabcs
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_allcl_chm                         ! Bc on all classes
  integer(ip),   pointer                :: &
       kfl_initi_chm(:),                   & ! Initial condition
       kfl_usrbc_chm(:)                      ! user boundary condition       
  real(rp),      pointer                :: &
       xinit_chm(:,:),                       & ! Initial pvalue parameter
       panat_chm(:,:)                        ! Natural bc parameters
  type(bc_nodes), pointer               :: &     
       tncod_chm(:)                          ! Node code type
  type(bc_nodes), pointer               :: &     
       tgcod_chm(:)                          ! Geometrical node code type
  type(bc_bound), pointer               :: &     
       tbcod_chm(:)                          ! Boundary code type
!--END REA GROUP
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  integer(ip),   pointer                :: &
       kfl_fixno_chm(:,:),                 & ! Nodal fixity 
       kfl_fixbo_chm(:,:)                    ! Boundary fixity
  real(rp),      pointer                :: &
       bvess_chm(:,:),                     & ! Essential bc values
       avtem_chm(:),                       & ! Time-averaged temperature
       avcon_chm(:),                       & ! Time-averaged concentration
       avvar_chm(:),                       & ! Time-averaged variance of RPV (VRPV)
       avime_chm(:),                       & ! Time-averaged normalized enthalpy
       avchm_chm(:),                       & ! Time-averaged chemical heat
       avmix_chm(:),                       & ! Time-averaged mixture fraction
       avmi2_chm(:),                       & ! Time-averaged squared of mixture fraction f*f
       avco2_chm(:)                          ! Time-averaged squared of concentration c*c

  integer(ip)                           :: &
       kfl_robin_chm,                      & ! Robin condition exists
       ncomp_chm,                          & ! Number of components 
       kfl_grdif_chm,                      & ! If there are gradients of conductivity
       kfl_tiaor_chm,                      & ! Original time accuracy
!       kfl_stead_chm,                      & ! Steady-state has been reached
       kfl_goite_chm,                      & ! Keep iterating
       iclas_chm,                          & ! Current class being solved
       iclai_chm,                          & ! Initial class
       iclaf_chm,                          & ! Final class
       ittot_chm,                          & ! Total number of iterations
       nskyl_chm,                          & ! Size of the skyline matrix
       kfl_goit2_chm,                      & ! Internal goite
       kfl_gocla_chm,                      & ! Intra-species iterations goite
       kfl_under_chm,                      & ! # undershoots
       kfl_overs_chm                         ! # overshoots

  real(rp)                              :: &
       dtinv_chm,                          & ! 1/dt
       dtcri_chm,                          & ! Critical time step
       resid_chm,                          & ! Residual for outer iterations
       pabdf_chm(10),                      & ! BDF factors
       rtpts_chm,                          & ! Global inner residual
       comin_chm,                          & ! Minimum concentration
       comax_chm,                          & ! Maximum concentration
       cputi_chm(10),                      & ! CPU time
       dtmat_chm,                          & ! Matrix time step
       xvael_chm(100),                     & ! Values
       grnor_chm                             ! gravity modulus 

  real(rp),     pointer                 :: &
       amatr_chm(:),                       & ! Matrix for ODE's
       rhsid_chm(:),                       & ! RHS for ODE's
       ripts_chm(:),                       & ! Class inner residual
       vmass_chm(:),                       & ! Mass matrix
       smatr_chm(:),                       & ! Constant matrix
       shsid_chm(:),                       & ! Constant RHS
       proje_chm(:,:)                        ! Projection
  integer(ip),  pointer                 :: &
       iarea_chm(:),                       & ! IA CSR format
       jarea_chm(:),                       & ! JA CSR format
       iskyl_chm(:),                       & ! DOE's skyline list
       idiag_chm(:),                       & ! DOE's skyline diagonal
       idima_chm(:)                          ! Position of diagonal in sparse matrix

  logical(lg)                            ::&
       kfl_rsta2_chm                         !restarted file was BDF
  
  !------------------------------------------------------------------------
  ! METEO model
  !------------------------------------------------------------------------

  real(rp)                              :: &
       tmete_chm,                          & ! Initial time for METEO model (properties)
       tsour_chm                             ! Initial time for METEO model (sources)
  real(rp),     pointer                 :: &
       veloc_chm(:,:),                     & ! Velocity          (from meteo file)
       densi_chm(:),                       & ! Density           (from meteo file)
       tempe_chm(:),                       & ! Temperature       (from meteo file)
       vfric_chm(:),                       & ! Friction velocity (from meteo file)
       hepbl_chm(:),                       & ! BL height         (from meteo file)
       walld_chm(:),                       & ! Wall distance     (from meteo file)
       lmoni_chm(:),                       & ! M-O length        (from meteo file)
       tmrat_chm(:,:),                     & ! Total Mass rate   (from meteo file)
       vterm_chm(:,:),                     & ! Terminal velocity per class
       accum_chm(:),                       & ! Ground Mass
       treac_chm(:,:),                     & ! Reactive term
       qvaporef_chm(:),                    & ! Qvapor           (from meteo file)
       qcloudref_chm(:),                   & ! Qvapor           (from meteo file)
       qrainref_chm(:)                       ! Qvapor           (from meteo file)
  
  !------------------------------------------------------------------------
  ! MECHANO-BIOLOGY model
  !------------------------------------------------------------------------

  real(rp),     pointer                 :: &
       proad_chm(:)                          ! Poteine adsorption
  !
  ! ADR type
  !
  type(ADR_typ)                         :: &
       ADR_read                              ! ADR type for reading

  type(ADR_typ), allocatable, target    :: & 
       ADR_chm(:)                            ! ADR type for nclas_chm

end module def_chemic
