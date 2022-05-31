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
  use mod_interp_tab,   only : typ_tab_coord, typ_lookup_table, typ_lookup_framework 
  use mod_ADR,          only : ADR_typ
#ifdef CANTERA
  use cantera
#endif

  implicit none
  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip), parameter ::                  &
       lun_times_chm   = 1911,               &  ! Time step info
       lun_time2_chm   = 1912,               &  ! Time step targets
       lun_resu1_chm   = 1921,               &  ! Result 1
       lun_resu2_chm   = 1922,               &  ! Result 2
       lun_resu3_chm   = 1923,               &  ! Result 3
       lun_resou_chm   = 1931,               &  ! Meteo source file
       lun_spcvg_chm   = 1932,               &  ! Species convergence file
       lun_droplet_chm = 1941                   ! Droplet results unit
       
  character(150)                        :: &
       fil_times_chm,                      &  ! Time step    
       fil_time2_chm,                      &  ! Time step target
       fil_resou_chm,                      &  ! Meteo
       fil_droplet_chm                        ! Droplet results file
 
  real(rp),      parameter :: &
       zepts = epsilon(1.0_rp)
  integer(ip),   parameter              :: &
       npara_chm=10,                       &   ! # parameters for sets
       npart_chm=10                            ! # temperature parameters
!--BEGIN REA GROUP
  !------------------------------------------------------------------------
  ! Physical problem: read in chm_reaphy
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_model_chm,                      & ! PDE-ODE model
       kfl_timei_chm,                      & ! Existence of du/dt
       kfl_advec_chm,                      & ! Existence of (a.grad)u
       kfl_diffu_chm,                      & ! Existence of -div[k. grad(u)]
       kfl_transport_chm,                  & ! Strategy for transport properties
       kfl_norma_chm,                      & ! Normalize mass fractions
       kfl_cotur_chm,                      & ! Flag to identify the RANS model for the Flamelet combustion model
       kfl_premix_chm,                     & ! Flag to identify the combustion model: premixed or non-premixed 
       kfl_varYc_chm,                      & ! Flag to define Yc variance model
       kfl_varZ_chm,                       & ! Flag to define Z variance model
       kfl_wallc_chm,                      & ! Flag to impose a zero source term at walls for the Flamelet combustion model
       kfl_radia_chm,                      & ! Flag to activate radiation model for Flamelet model
       kfl_pfa_chm,                        & ! Flag to activate PFA reduction model
       kfl_key_chm,                        & ! Flag to determine number of key species 
       kfl_cont_chm,                       & ! Number of controlling varibale TDAC
       kfl_z_chm,                          & ! Flag to determine z computation type
       kfl_freq_chm,                       & ! Number of timesteps between DAC reductions 
       kfl_tdac_write_chm,                 & ! Write Reduced table 
       kfl_spray_chm,                      & ! Flag to activate spray model
       kfl_ufpv_chm,                       & ! Unsteady Flamelet Progress Variable model
       kfl_heat_loss_chm,                  & ! Heat loss model
       kfl_lookg_chm,                      & ! Lookup on Gauss points 
       kfl_tab_fw_chm,                     & ! Index of table lookup framework
       kfl_post_fw_chm,                    & ! Index of postprocessing table lookup framework
       kfl_spec_name_chm,                  & ! Number of Fields (Finite Rate)
       kfl_entropy_chm,                    & ! Activate entropy stable viscosity stabilization method
       kfl_field_chm(3),                   & ! Flag to activate the initialization by fields (1) Species fields, (2) Temperature field, (3) alpha (CMC model)
       nclas_chm,                          & ! Number of particle classes
       nspec_chm                             ! Total number of species

  !
  ! Flamelet model variables
  ! 
  type(typ_tab_coord), pointer  :: &
         chi_shape_coords(:),  & ! control variable discretization for shape function of scalar dissipation rate
         yc_scaling_coords(:), & ! control variable discretization for progress variable limits 
         h_scaling_coords(:),  & ! control variable discretization for enthalpy limits
         posttab_coords(:),    & ! control variable discretization for post processing lookup table
         table_coords(:)         ! control variable discretization for lookup table

  type(typ_lookup_table), pointer :: & 
      chi_shape_tab,             & ! Table for shape function of scalar dissipation rate
      yc_scaling_tab,            & ! Table for progress variable limits
      h_scaling_tab,             & ! Table for enthalpy limits
      posttab_tab,               & ! Table for post processing lookup table
      table_tab                    ! Table for lookup table 

  type(typ_lookup_framework), pointer :: &
      table_fw,                  & ! Framework for lookup table
      posttable_fw                 ! Framework for lookup table


  ! Flamelet properties on Gauss points
  type(r3p),pointer        :: &
       mass_gp(:),  &             ! mass source term 
       hk_gp(:),    &             ! enthalpy at gauss points (sensible + chemical)
       rrt_gp(:)                  ! subscale source for variance

  !
  ! Finite-rate combustion model
  !
  real(rp),  pointer                    :: & 
       grad_Yk(:,:,:),                     & ! Gradient of species Yk
       grad_Yk_hk(:,:,:),                  & ! Gradiente of species Yk times the enthalpy (chem.+sens.) of species k
       grad_T(:,:),                        & ! Gradient of temperature
       grad_h(:,:),                        & ! Gradient of enthalpy (chem.+sens.)
       aux_nodes(:),                       & ! Auxiliary variable at nodes
       coeff_cp_k(:,:),                    & ! NASA polinomial coefficients for individual species
       W_k(:),                             & ! Species molecular weights
       Y_k_n(:,:,:),                       & ! Local mass fractions of all species for PFA reduction methods
       enthalpy_transport_nodes(:,:),      & ! Enthalpy transport by diffusion
       grad_k_single(:,:),                 & ! Gradient related to a given field (e.g. species, etc.). It is used as an auxiliary variable.
       Le_k(:),                            & ! Lewis number for each species
       hrr_chm(:),                         & ! Heat release field (W/m3)
       hrr_avg_chm(:)                        ! Averaged heat release field (W/m3)
  
  !
  ! CMC turbulent combustion model
  !

  integer(ip), pointer                  :: &
      kfl_bc_type_spec_CMC_chm(:)            ! It values 1 if for all the physical points the boundary conditions for all the chemical species are of the same type and 0 otherwise


  real(rp),  pointer                    :: & ! The order of the dimensions is (mixt fraction, points, variables)
      enthalp_CMC_chm(:,:),                & ! Conditional enthalpy for CMC h at nodes (nZ,npoin)
      temp_CMC_chm(:,:),                   & ! Conditional temperature for CMC at nodes (nZ,npoin)
      Yk_CMC_chm(:,:,:),                   & ! Conditional mass fractions for CMC Y_k at nodes (nZ,npoin,nclas)
      src_Yk_CMC_chm(:,:,:),               & ! Conditional chemical source term (nZ,npoin,nclas)
      Yk_int_CMC_chm(:,:),                 & ! Unconditional mass fractions for CMC Y_k at nodes (npoin,nclas)
      densi_int_CMC_chm(:),                & ! Density computed from CMC to be transferred to CFD (npoin)
      visco_lam_int_CMC_chm(:),            & ! Laminar viscosity from CMC to be transferred to CFD (npoin)
      enthalp_int_CMC_chm(:),              & ! Unconditional enthalpy for CMC Y_k at nodes (npoin)
      temp_int_CMC_chm(:),                 & ! Unconditional temperature for CMC Y_k at nodes (npoin)
      src_Yk_int_CMC_chm(:,:),             & ! Unconditional chemical source term (npoin,nclas)
      hrr_CMC_chm(:,:),                    & ! Heat release (nZ,npoin)
      veloc_CFD_chm(:,:),                  & ! Velocity field from CFD calculation (ndime,npoin)
      Zavg_CFD_chm(:),                     & ! Mixture fraction field from CFD (npoin)
      Zvar_CFD_chm(:),                     & ! Mixture fraction variance field from CFD (npoin)
      Xtot_CFD_chm(:),                     & ! Total scalar dissipation rate from CFD (npoin) (sgs+solved)
      grad_Zavg_CFD_chm(:,:),              & ! Gradient of mixture fraction field from CFD (npoin,ndime)
      visco_turb_CFD_chm(:),               & ! Mass turbulent diffusion coefficient from CFD (npoin)
      deriv2_Yk_CMC_chm(:,:,:),            & ! Second derivative in mixture fraction direction for the conditional species (nZ,npoin,nclas)
      deriv2_enthalp_CMC_chm(:,:),         & ! Second derivative in mixture fraction direction for the conditional enthalpy (nZ,npoin)
      bvess_CMC_chm(:,:,:),                & ! Essential bc conditional values for CMC (nZ,npoin,nvar_CMC_chm)
      bvess_ufield_CMC_chm(:,:),           & ! Essential bc unconditional fields (npoin,nvar_CMC_chm)
      T_bc_CMC_chm(:),                     & ! Temperature at boundaries if kfl_weigh_in_eq_CMC_chm activated
      react_scal_bc_CMC_chm(:,:),          & ! Reactive scalars (enthalpy and species) at the Z=0 and Zs if kfl_weigh_in_eq_CMC_chm activated
      Z_CMC_chm(:),                        & ! Mixture fraction vector for CMC (nZ)
      diff_Z_CMC_chm(:),                   & ! Vector with the increments in mixture fractions for Z_CMC_chm (nZ-1)
      Xnormalized_prof_CMC_chm(:),         & ! Scalar dissipation rate normalized profile (nZ)
      rscal_inert_CMC_chm(:,:),            & ! Inert profile for reactive scalars (enthalpy and species) (nZ,nclas+1)
      rscal_equil_CMC_chm(:,:),            & ! Equilibrium profile for reactive scalars (enthalpy and species) (nZ,nclas+1)
      temp_inert_CMC_chm(:),               & ! Inert profile for temperature (nZ)
      temp_equil_CMC_chm(:),               & ! Equilibrium profile for temperature (nZ)
      Z_AMC_CMC_chm(:),                    & ! Mixture fraction vector for AMC model (nZ_AMC_CMC_chm)
      S_AMC_CMC_chm(:),                    & ! Segregation factor vector for AMC model (nS_AMC_CMC_chm)
      Xintegrated_table_AMC_CMC_chm(:,:)     ! Table that contains the integrated profiles for scalar dissipation rate for AMC model from CMC model (nZ_AMC_CMC_chm, nS_AMC_CMC_chm)


  integer(ip)                           :: &
      kfl_weigh_in_eq_CMC_chm,             & ! When starting from scratch activate the option to find the initial solution from weighing the inert and equilibrium solutions
      kfl_split_CFD_CMC,                   & ! Flag to split CFD and CMC into two different executions: 0 CFD and CMC in same execution and 1 if CFD and CMC separated (MODIFY: THIS SHOULD BE IN THE KERNEL). If 1 make turbulent diffusion matrix points to nu_t
      kfl_solve_enth_CMC_chm,              & ! 0 if enthalpy is not solved and 1 if it is transported
      kfl_start_CMC_chm = 1,               & ! 1 for the initial time step; 0 otherwise
      nZ_CMC_chm,                          & ! Number of slices in mixture fraction space for CMC conditioning. Due to different reasons when using finite rate we take nZ_CMC_chm = 3
      nZ_AMC_CMC_chm,                      & ! Number of mixture fractions for scalar dissipation rate integration (AMC model)
      nS_AMC_CMC_chm,                      & ! Number of segregation factors for scalar dissipation rate integration (AMC model)
      nsize_mf_vec_chm,                    & ! Length of the mixture fraction path
      imixf_rk,                            & ! Mixture fraction iterator in RK scheme
      nvar_CMC_chm,                        & ! Number of variables to be solved in CMC: nvar_CMC_chm = nclas_chm+1 (species + enthalpy)
      index_N2                               ! Index in the mechanism for N2


  real(rp)                              :: &
      Zs_CMC_chm,                          & ! Saturation mixture fraction
      Smax_AMC_CMC_chm,                    & ! Maximum segregation factor for AMC model
      S_threshold                            ! Threshold segregation factor for integrations


  type(r3p),pointer                     :: &
       ! Following variables save conditional values
       condu_gp_CMC_chm(:),                & ! Thermal conductivity at Gauss points (nelem)%(pgaus,1,1)
       sphec_gp_CMC_chm(:),                & ! Specific heat at Gauss points (nelem)%(pgaus,1,1)
       visco_gp_CMC_chm(:),                & ! Viscosity at Gauss points (nelem)%(nZ,pgaus,1)
       spvol_gp_CMC_chm(:)                   ! Specific volume at Gauss points (nelem)%(nZ,pgaus,1)


  character(len=:), allocatable         :: &
       mf_vec_path_CMC_chm                   ! Mixture fraction vector path for CMC

  !
  ! End of CMC turbulent combustion model variables
  !


  character(len=:), allocatable ::         &
       mechanism_path,                     & ! Mechanism name for Finite Rate chemistry
       Red_spec                              ! Species to be reduced 

  integer(ip)                           :: &
       nsize_mech_name,                    & ! Number of characters mechanism name
       nsize_red                             ! Number of characters of reduced Species 
  !
  ! Global variables
  !
  real(rp),  pointer                    :: & 
       rspec_chm(:,:),                     & ! Flamelet model: species mass fraction for radiation model
       zgradmax_chm(:),                    & ! Maximum gradient of the mixture fraction that defines phi_chm=1
       phi_chm(:),                         & ! Weighting factor for the hybrid model
       dt_rho_chm(:),                      & ! Projection of dt/rho
       dt_chm(:),                          & ! Projection of dt for interface tracking in ELSA model
       avden_chm(:),                       & ! Time-averaged density
       avY_chm(:),                         & ! Time-averaged reaction progress variable Yc or C
       avY2_chm(:),                        & ! Time-averaged reaction progress variable squared Yc*Yc or C*C
       avZv_chm(:),                        & ! Time-averaged variance of Z
       avchm_chm(:),                       & ! Time-averaged chemical heat
       avZ_chm(:),                         & ! Time-averaged mixture fraction Z
       avYv_chm(:),                        & ! Time-averaged variance reaction progress Yc 
       avmsk_chm(:),                       & ! Time-averaged mass source term from spray 
       xYr_chm(:),                         & ! Scalar dissipation rate of Yc (resolved part)
       xZr_chm(:),                         & ! Scalar dissipation rate of Z  (resolved part)
       xYs_chm(:),                         & ! Scalar dissipation rate of Yc (subgrid part)
       xZs_chm(:),                         & ! Scalar dissipation rate of Z  (subgrid part)
       avxYr_chm(:),                       & ! Average scalar dissipation rate of Yc (resolved part)
       avxZr_chm(:),                       & ! Average scalar dissipation rate of Z  (resolved part)
       avxYs_chm(:),                       & ! Average scalar dissipation rate of Yc (subgrid part)
       avxZs_chm(:),                       & ! Average scalar dissipation rate of Z  (subgrid part)
       avZ2_chm(:),                        & ! Time-averaged squared of mixture fraction Z*Z
       avposttab_chm(:,:)                    ! Time-averaged tabulated postprocessing quantities

  !
  ! Spray variables
  !
  real(rp),  pointer                    :: & 
       avL_chm(:),                         & ! Time-averaged liquid volume fraction phi_L
       avL2_chm(:),                        & ! Time-averaged liquid volume fraction squared phi_L*phi_L
       Sigma_chm(:),                       & ! Interface surface density Sigma
       Sigm0_chm(:),                       & ! Interface surface density Sigma_0 or Sigma_min
       d32_chm(:),                         & ! Sauter mean diameter
       avS_chm(:),                         & ! Time-averaged interface surface density Sigma
       avS0_chm(:),                        & ! Time-averaged interface surface density Sigma_0 or Sigma_min
       avd32_chm(:)                          ! Time-averaged Sauter mean diameter

  type(r3p),  pointer                   :: &
       sigma_gp_chm(:),                    & ! Interface surface density Sigma
       d32_gp_chm(:),                      & ! Sauter mean diameter
       dummy_enthalpy(:),                  & ! enthalpy at gauss points (sensible + chemical)
       sigma0_gp_chm(:)                      ! Interface surface density Sigma_0 or Sigma_min
  !
  ! Level set variables for spray
  !
  real(rp),  pointer                    :: &
       lap_phi_levSet_chm(:),              & ! Laplacian of the gradient of level set function
       grad_phic_levSet_chm(:,:),          & ! Gradient of phi*(1-phi)
       grad_phi_levSet_chm(:,:)              ! Gradient of level set function phi

  !
  ! Sectional soot model
  !
  integer(ip)                           :: &
       kfl_soot_chm                        ! Activation soot model
!!DMM       indexS(5),                          & ! Index of soot sources models 
!!DMM       gasCoupling_ssm,                    & ! Coupling with gas phase
!!DMM       nPAH,                               & ! Number of species involved in nucleation = condensation 
!!DMM       nSurf,                              & ! Number of species involved in surface growth
!!DMM       idNucl(10),                         & ! Index species involved in nucleation 
!!DMM       idCond(10),                         & ! Index species involved in condensation
!!DMM       idSurf(10),                         & ! Index species involved in surface growth
!!DMM       ID_NUCL,                            & ! Index nucleation process
!!DMM       ID_COND,                            & ! Index condensation process
!!DMM       ID_COAG,                            & ! Index coagulation process
!!DMM       ID_SURF,                            & ! Index surface growth process
!!DMM       nsect_chm                             ! Number of sections for soot sectional model

  !
  ! Others
  !       
  real(rp)                              :: &
       radwt_chm,                          & ! Wall temprature for radiation model
       dac_crit_chm,                       & ! Flag for crit value in PFA
       dac_cor_chm,                        & ! Corrleation thershold for CODAC
       bf_fuel_chm,                        & ! Bilgers Formula Fuel 
       bo_oxy_chm,                         & ! Bilgers Formula Oxy 
       sorad_chm,                          & ! Source radius
       socen_chm(3),                       & ! Source center
       prthe_chm,                          & ! Thermodynamic pressure (if NASTIN not activated)
       surf_tension_chm,                   & ! Surface tension for sprays
       hrr_int_chm                           ! Total heat release in the domain (W)

  integer(ip)                           :: &
       nreac_chm,                          & ! Number of reactions
       stofu_chm(3)                          ! Which is fuel and which is oxygen for equivalence ratio

  !
  ! Variables transmitted to nodes later
  !
  integer(ip), pointer                  :: &
       React_ind(:,:),                     & ! Local reduced mechanims based on PFA
       Field_ind_chm(:),                   & ! Index of field variable 
       Corr_chm(:)                           ! Local correlation for dynamic reduction 

  real(rp),  pointer                    :: &      
       diffu_chm(:,:)                        ! Diffusion constants

  real(rp),  pointer                    :: & 
       entha_chm(:,:),                     & ! Enthalpy (sensible + chemical) for each species
       elem_h(:),                          & ! Elemental mass fraction of H 
       elem_c(:),                          & ! Elemental mass fraction of C
       elem_o(:),                          & ! Elemental mass fraction of O
       elem_n(:),                          & ! Elemental mass fraction of N
       mixfr_chm(:),                       & ! Mixture Fraction Chemic (finite rate)
       prog_var_chm(:),                    & ! Progress Variable TDAC
       sum_reac_chm(:),                    & ! Sum of Reactions at a node 
       src_chm(:,:)                          ! Chemical Source - Finite Rate

  integer,parameter                     :: &
       maxsp_chm = 1000                       ! Max number of species tracked in combustion code

  !
  ! Eulerian droplet identification variables
  !
  integer(ip)                           :: &
       kfl_droplet_id_chm,                 & ! Droplet identification flag           
       ndrop_chm,                          & ! Total number of identified Eulerian droplets
       droplet_postprocess_frequency_chm     ! Droplet postprocess frequency

  real(rp)                              :: &
       levelSet_threshold_chm,             & ! Level Set threshold for droplet identification
       droplet_compactness_limit_chm,      & ! Compactness value below which a cluster won't be considered as a droplet
       droplet_max_diameter_chm,           & ! Max. diameter above which a cluster won't be considered as a droplet
       droplet_h_factor_chm                  ! Mesh size factor to define a max. diameter

  real(rp),  pointer                    :: &
       volume_drop_chm(:),                 & ! Droplet volume
       diameter_drop_chm(:),               & ! Droplet diameter
       centroid_drop_chm(:,:),             & ! Droplet centroid
       compactness2_drop_chm(:),           & ! Droplet compactness_2
       volume_cluster_chm(:)                 ! Cluster volume

  !------------------------------------------------------------------------
  ! Numerical problem: read in chm_reanut
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_ellen_chm,                      & ! =0,1 for min/max element length
       kfl_taust_chm,                      & ! Tau calculation option
       kfl_shock_chm,                      & ! Shock capturing type 
       kfl_stabi_chm,                      & ! Stabilization strategy
       kfl_limit_chm,                      & ! Limiter
       kfl_tiacc_chm,                      & ! Temporal accuracy
       kfl_tibub_chm,                      & ! Time integration of bubble
       kfl_split_chm,                      & ! Splitting algorithm
       kfl_tisch_chm,                      & ! Time integration scheme
       kfl_normc_chm,                      & ! Norm of convergence
       kfl_dtcri_chm,                      & ! dt criteria
       kfl_negat_chm,                      & ! Startegy for negative concentrations
       kfl_posit_chm,                      & ! Startegy for too positive concentrations
       kfl_warni_chm,                      & ! Warn about points with zero sum mass
       kfl_temli_chm                         ! Flag to activate a T limiter to compute reaction rates

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
       relax_chm                             ! Relaxation of update

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
       kfl_conbc_chm,                      & ! Constant boundary conditions
       kfl_allcl_chm,                      & ! Bc on all classes
       kfl_fields_scale_chm                  ! Initialization by unscale (== 1) or scale (== 0) values
  integer(ip),   pointer                :: &
       kfl_initi_chm(:)                      ! Initial condition
  real(rp),      pointer                :: &
       xinit_chm(:,:)                        ! Initial pvalue parameter
  type(bc_nodes), pointer               :: &     
       tncod_chm(:)                          ! Node code type
  type(bc_nodes), pointer               :: &     
       tgcod_chm(:)                          ! Geometrical node code type
  type(bc_bound), pointer               :: &     
       tbcod_chm(:)                          ! Boundary code type
!--END REA GROUP

  !
  ! Others
  !
  integer(ip),   pointer                :: &
       kfl_fixno_chm(:,:),                 & ! Nodal fixity 
       kfl_fixbo_chm(:,:),                 & ! Boundary fixity
       kfl_funno_chm(:,:),                 & ! Function # for node BC
       kfl_funtn_chm(:,:)                    ! Function type for node BC

  real(rp),      pointer                :: &
       bvess_chm(:,:)                      ! Essential bc values

  integer(ip)                           :: &
       kfl_robin_chm,                      & ! Robin condition exists
       ncomp_chm,                          & ! Number of components 
       kfl_grdif_chm,                      & ! If there are gradients of conductivity
       kfl_tiaor_chm,                      & ! Original time accuracy
       kfl_goite_chm,                      & ! Keep iterating
       iclas_chm,                          & ! Current class being solved
       iclai_chm,                          & ! Initial class
       iclaf_chm,                          & ! Final class
       ittot_chm,                          & ! Total number of iterations
       nskyl_chm,                          & ! Size of the skyline matrix
       kfl_goit2_chm,                      & ! Internal goite
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
       shsid_chm(:)                          ! Constant RHS

  integer(ip),  pointer                 :: &
       iarea_chm(:),                       & ! IA CSR format
       jarea_chm(:),                       & ! JA CSR format
       iskyl_chm(:),                       & ! DOE's skyline list
       idiag_chm(:),                       & ! DOE's skyline diagonal
       idima_chm(:)                          ! Position of diagonal in sparse matrix

  logical(lg)                            ::&
       kfl_rsta2_chm                         ! Restarted file was BDF
  
  type(ADR_typ), allocatable, target    :: & 
       ADR_chm(:)                            ! ADR type for nclas_chm

#ifdef CANTERA
  type(phase_t) gas_chm
#endif


end module def_chemic
