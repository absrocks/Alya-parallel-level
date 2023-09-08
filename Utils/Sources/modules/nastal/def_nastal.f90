!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    def_nastal.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Definition module
!> @details Definition module
!> @} 
!-----------------------------------------------------------------------
module def_nastal
!-----------------------------------------------------------------------
!    
! Heading for the module subroutine
!
!-----------------------------------------------------------------------
  use def_kintyp

!------------------------------------------------------------------------
! Parameters
!------------------------------------------------------------------------

  integer(ip), parameter ::&
       lun_chkpo_nsa = 605, lun_pro2d_nsa = 606,&
       lun_force_nsa = 607, lun_maxmi_nsa = 610,&
       lun_dumb1_nsa = 620, lun_dumb2_nsa = 621
  real(rp), parameter :: &
       zensa = epsilon(0.0_rp)
  integer(ip), parameter ::&                     ! # postprocess variables
       nvecp_nsa=20,       &
       nspep_nsa=50,       &
       nscap_nsa=40,       &
       ncoef_nsa=10                              ! # coefficient for properties
  integer(ip), parameter ::&                     ! # max sets
       msets_nsa=10

!------------------------------------------------------------------------
! Physical problem: read, modified or defined in nsa_reaphy
!------------------------------------------------------------------------
  ! nsa_reaphy types

  type geome_posit_nsa
     integer (ip) :: kflag= 0
     integer (ip) :: geometry = 0
     real(rp)     :: value= 0.0_rp
     real(rp)     :: center(3)= 0.0_rp
     real(rp)     :: radius = 0.0_rp
  end type geome_posit_nsa

  type(geome_posit_nsa), save :: iposi_nsa(8)

!--BEGIN REA GROUP
  integer(ip) ::&
       kfl_timei_nsa,&              ! Existence of du/dt
       kfl_stead_nsa,&              ! Steady-state has been reached 
       kfl_advec_nsa,&              ! Existence of (u.grad)u
       kfl_visco_nsa,&              ! Viscous term 
       kfl_inico_nsa,&              ! Initial conditions type
       kfl_infun_nsa,&              ! Special initial fields functions
       kfl_benme_nsa,&              ! Particular meteo benchmark problem
       kfl_spher_nsa,&              ! pseud-3D cylinder/sphere flag
       kfl_physics_nsa,&            ! Flag to turn on and off the microphysics
       kfl_botta_nsa,&              ! Botta anfd Klein equilibrium case
       kfl_rearst_nsa,&             ! Flag to read local restart file
       kfl_sptyp_nsa,&              ! Sponge type (Lilly and Klemp 1978 or simpler Giraldo's)
       kfl_hysta_nsa,&              ! Hydrostatic correction
       kfl_brunt_nsa,&              ! Brunt frequency flag
       kfl_inifi_nsa(3),&           ! Initial fields
       kfl_inibu_nsa   ,&           ! Built-in initial fields
       kfl_inkee_nsa(5),&           ! What kind of initial fields are used. Order: velocity, density, temperature, pressure and Brunt
       kfl_inleb_nsa,&
       kfl_dampf_nsa,&              ! Van Driest near-wall damping function for Smagorisnky model
       kfl_mfrco_nsa,&              ! Mass flow rate control to keep a constant bulk velocity (Ub) on a given set
       kfl_refpdt_nsa,&             ! Reference p-d-t: compute pressure from t and rho as default
       kfl_iposi_nsa,&              ! Positional initial fields
       kfl_cylind_nsa = 1           ! By default the thermal perturbation is a circle/ellipses in 2D, or a cylinder in 3D. 

  integer(ip) ::         &
       nmate_nsa,      &         ! # of materials
       ndofn_nsa,      &         ! # of d.o.f. of the problem
       ndof2_nsa,      &         ! ndofn_*ndofn_
       nkeep_nsa,      &         ! ndofn_nsa + 2, size of the base state for KEEPHYDRO
       nevat_nsa,&               ! Element matrix dimension = ndofn_nsa*mnode
       nevab_nsa,&               ! Element velocity dof=ndime*nnode
       nflub_nsa,&               ! Number of equations where flux boundary conditions could be imposed
       ncomp_nsa,&               ! Number of components of the fields
       nfrap_nsa,&               ! Number of fractional steps for FRK
       nromp_nsa,&               ! Number of components of the residuals
       ntomp_nsa,&               ! Number of components of the temporal internal variable
       nfiel_nsa(10),&           ! Fields assignement
       mfreq_nsa,&               ! Base number of frequency mode for subscales
       frmax_nsa                 ! Maximum number of frequency mode for subscales

  real(rp) ,     pointer      :: &
       frequ_nsa(:)              ! Nodal frequency field of Nondiagonal VMS

  real(rp) ,     pointer      :: &
       crens_nsa(:,:)              ! Convergence field

  integer(ip) ::&
       ivert_nsa,&                  ! Vertical coordinate from gravity definition
       lawde_nsa,&                  ! Law for rho
       lawvi_nsa,&                  ! Law for mu
       lawpo_nsa,&                  ! Law for sig    
       lawst_nsa,&                  ! State law    
       nuinlet_nsa,&                ! Nodes with total pressure condition
       mfrse_nsa                    ! Set from which the mass flow rate is based on

  real(rp) ::&
       dtinv_nsa,&                  ! 1/dt
       adgam_nsa,&                  ! Adiabatic exponent
       cpcoe_nsa,&                  ! Cp coefficient
       cvcoe_nsa,&                  ! Cv coefficient
       rgasc_nsa,&                  ! R gas constant
       rgacv_nsa,&                  ! R / Cv
       vispa_nsa(ncoef_nsa),&       ! Parameters for the law of viscosity
       fcons_nsa    ,&              ! Convection term
       fvins_nsa    ,&              ! Viscous term
       entro_nsa,&                  ! Reference inflow entropy
       tstag_nsa,&                  ! Stagnation temperature
       trefa_nsa,&                  ! Temperature recovery factor
       cppra_nsa,&                  ! Cp / Prandtl
       cpprt_nsa,&                  ! Cp / Turbulent Prandtl
       pbaro_nsa,&                  ! Reference pressure at ground level for barometric law
       xbubb_nsa(6,5),&             ! Central position for 5 bubbles and their radii
       teano_nsa(5)  ,&             ! Temperature anomaly (amplitude)
       envol_nsa     ,&             ! Energy integral over the total volume
       devol_nsa     ,&             ! Density integral over the total volume (i.e. total mass)
       sofac_nsa     ,&             ! Sound speed factor
       turbu_nsa                    ! Eddy viscosity coefficient

!
! Physical properties and reference values for coupled problems (CHEMIC)
! 1. moist air (total)
! 2. water vapor
! 3. cloud water
! 4. liquid water
!
  real(rp) ::&
       mixrt_nsa,   &               ! Total mixing ratio (constant for now)
       rgcou_nsa(5),&               ! R gas constant for species (COUPLING WITH CHEMIC)
       cpcou_nsa(5),&               ! Cp constant for species (COUPLING WITH CHEMIC)
       cvcou_nsa(5),&               ! Cv constant for species (COUPLING WITH CHEMIC)
       lacou_nsa(5),&               ! Reference latent heats (COUPLING WITH CHEMIC)
       pacou_nsa(10)                ! Coupling parameters (COUPLING WITH CHEMIC)
!
! Physical properties and reference values
!
  real(rp) ::&
       llrho_nsa(4),&           ! Density for lax-liu
       llvel_nsa(3,4),&         ! Velocity for lax-liu
       llpre_nsa(4),&           ! Density for lax-liu
       llcen_nsa(3),&           ! Center for the lax-liu
       gravi_nsa(3) ,&          ! Gravity vector
       grnor_nsa    ,&          ! Gravity norm
       veloc_nsa(3) ,&          ! Reference velocity (rho)
       densi_nsa,&              ! Reference density (rho)
       tempe_nsa,&              ! Reference temperature (rho)
       mixrv_nsa,&              ! Dummy value for the vapor mixing ratio (never really used for any def. or computation).
       speed_nsa,&              ! Reference velocity module (u)
       press_nsa,&              ! Reference pressure (p_0)
       press_dynamic_nsa,&      ! Reference dynamic pressure (q_0)
       visco_nsa,&              ! Reference dynamic viscosity (mu)
       gravm_nsa(5,5),&         ! Gravity matrix
       mfrgr_nsa,&              ! Mass flow rate growth parameter
       mfrub_nsa,&              ! Bulk velocity Ub target
       ubpre_nsa,&              ! Bulk velocity Ub target from previous time-step
       thdif_nsa,&              ! Reference thermal diffusion (k)
       prand_nsa,&              ! Prandtl number (SUPERSEDES k)
       pratu_nsa,&              ! Turbulent Prandtl number
       rreyn_nsa,&              ! Reference Reynolds number (SUPERSEDES mu)
       rmach_nsa,&              ! Reference Mach number (SUPERSEDES Cp)
       xsoun_nsa,&              ! Reference Sound speed
       poros_nsa(ncoef_nsa),&   ! Porosity (sig)
       cleng_nsa,&              ! Reference characteristic length
       uinlet_nsa,&             ! 
       spein_nsa,&              ! 
       axyin_nsa,&              ! 
       axzin_nsa,&              ! 
       ayzin_nsa,&              ! 
       attyz_nsa,&              ! 
       attxy_nsa,&              ! 
       attxz_nsa,&              ! 
       brure_nsa,&              ! Brunt-vaissala frequency reference value
       mowei_nsa,&              ! Molecular weight of the mixture
       runiv_nsa                ! Universal gas constant Ro [J/K/mol]

  real(rp),     pointer     ::&
       conce_nsa(:,:,:)         ! SM Local concentration array conce.Allocated and used only when coupling occurs.

  integer(ip),  pointer     ::&
       brunt_nsa(:),&               ! Brunt-vaissala frequency (nodes)
       lmate_nsa(:),&               ! Materials (elements)
       lmatn_nsa(:)                 ! Materials (nodes)
  !
  !      Coupling
  !
  integer(ip) ::         &
       kfl_coupl_nsa,&                ! Coupling type
       kfl_cotur_nsa                  ! Coupling with turbulence (if <0 (Boussinesq approx) & if >0 (coupling with TURBUL) )   

  !------------------------------------------------------------------------
  ! Additional domain variables used only locally with Kessler (SM):
  !------------------------------------------------------------------------
  integer(ip) ::&
       kfl_adiff_nsa, &        ! Flag for artificial diffusion
       vtkrestart_nsa,&        ! Restart interval for Kessler VTK output
       kfl_ansou_nsa,&         ! Flag to discriminate, in initial conditions, between analytic or read sounding
       istep_nsa, &
       nvar_nsa,&              ! number of variables (dynamics+tracers)
       nelx_nsa,&              ! number of elements in x
       nely_nsa,&              ! number of elements in y !used only if ndime == 3, otherwise z is used in the vertical.
       nelz_nsa,&              ! number of elements in z
       nx_nsa,  &              ! number of nodes in x
       ny_nsa,  &              ! number of nodes in y !used only if ndime == 3, otherwise z is used in the vertical.
       nz_nsa,  &              ! number of nodes in z
       ncol_nsa,&              ! number of columns (used for Kessler).
       kfl_uniformvelocity_nsa,&
       kfl_thetac_nsa,&
       kfl_tracerc_nsa,&
       kfl_rc_nsa,&
       kfl_xr_nsa,&
       kfl_yr_nsa,&       
       kfl_zr_nsa,&
       kfl_xc_nsa,&
       kfl_yc_nsa,&
       kfl_zc_nsa
  
  real(rp)  ::   &             !Sponge parameters
        dxs_nsa, &
        dys_nsa, &
        dzs_nsa, &
        ampx_nsa,&
        ampy_nsa,&
        ampz_nsa,&
        xmin_nsa,&
        xmax_nsa,&
        ymin_nsa,&
        ymax_nsa,&
        zmin_nsa,&
        zmax_nsa,&
        xradi_nsa,&
        yradi_nsa,&
        zradi_nsa,&
        rc_nsa,   &
        xc_nsa,   &
        yc_nsa,   &
        zc_nsa,   &
        uvelo_nsa,&
        vvelo_nsa,&
        wvelo_nsa,&
        thetac_nsa,&
        tracerc_nsa

  integer(ip),  pointer      ::& 
       xcol_nsa(:),            &
       ycol_nsa(:),            &
       node_column_nsa(:,:),   &
       intma_column_nsa(:,:)

  real(rp)  ::&
       kdiff_nsa               ! Coefficient of artificial diffusion
  
  !------------------------------------------------------------------------
  ! Numerical problem: read, modified or defined in nsa_reanut
  !------------------------------------------------------------------------

  integer(ip) ::&
       kfl_timet_nsa,&              ! Implicit or explicit scheme
       kfl_adres_nsa,&              ! Adaptive subiterations tolerances
       kfl_matri_nsa,&              ! Explicit scheme but computing the rhs as matrix * vector
       kfl_algor_nsa,&              ! Monolithic algorithm flag
       kfl_diagi_nsa,&              ! Diagonal implicit terms active
       kfl_lotim_nsa,&              ! Diagonal local time steps terms active
       kfl_foreg_nsa,&              ! Forced regime: compressible or incompressible
       kfl_relat_nsa,&              ! Relativistic flows
       kfl_isent_nsa,&              ! Isentropic flows
       kfl_pertu_nsa,&              ! Using perturbation variables (meto only)
       kfl_theta_nsa,&              ! Using perturbation variables (meto only)
       kfl_ncons_nsa,&              ! Flag to use the Non-conservative set
       kfl_goite_nsa,&              ! Keep iterating
       kfl_stabi_nsa,&              ! Stabilization method
       kfl_stafl_nsa,&              ! Stabilization term formulation: using fluxes or using the Jacobian matrixes
       kfl_galer_nsa,&              ! Galerkin terms form: with jacobian A or with fluxes
       kfl_repro_nsa,&              ! Stabilization based on residual projection
       kfl_taufa_nsa(10,2),&        ! Tau and time stabilization factors
       kfl_shock_nsa(5),&           ! Shock capturing type 
       kfl_hconv_nsa,&              ! Speed vector correction for hconv
       kfl_penal_nsa,&              ! Penalization
       kfl_turbu_nsa,&              ! Turbulent flow
       kfl_weigh_nsa,&              ! Weighting of du/dt
       kfl_tiacc_nsa,&              ! Temporal accuracy
       kfl_taudi_nsa,&              ! Diagonal tau
       kfl_track_nsa,&              ! Subscales tracking
       kfl_higha_nsa,&              ! High aspect ratio elements (use a special strategy for computing tau and dt)
       kfl_resmo_nsa,&              ! Residual smoothing
       kfl_skews_nsa,&
       kfl_normc_nsa,&              ! Norm of convergence
       kfl_linea_nsa,&              ! Linearization (RHS=0, Picard=1, Newton=2)
       kfl_ximpl_nsa(4),&           ! Explicit terms in the implicit (Viscous, Convection)
       kfl_delun_nsa,&              ! Delta form
       kfl_algso_nsa,&              ! Type of algebraic solver
       kfl_tisch_nsa,&              ! Time integration scheme
       kfl_tiext_nsa,&              ! Externally fixed time step flag
       kfl_cfllo_nsa,&              ! Local CFL
       kfl_dttyp_nsa(10),&          ! Time increment type (global or local) per equation
       kfl_dtadj_nsa,&              ! Time increment type (global or local) adjustment
       kfl_fasts_nsa,&              ! Fast stationary (number of neighbors layers)       
       kfl_unkse_nsa,&              ! Conservative or primitive unknowns set
       kfl_rayle_nsa,&              ! Rayleigh damping
       kfl_lopre_nsa,&              ! Local preconditioning (if ACTIVE, kfl_unkse_nsa is forced to 1)
       kfl_zero_initial_velocity_nsa,&  ! Identifies those problems with zero initial velocity
       kfl_pseud_nsa,&              ! Pseudo time-step is used or not
       kfl_locti_nsa,&              ! Local time step 
       kfl_nacdr_nsa,&              ! Nastal cdr
       kfl_reate_nsa,&              ! Reate in tau and dt
       kfl_modfi_nsa,&              ! Modify kfl_fixno_nsa
       kfl_mod_elmop_nsa,&          ! Wether using elmoperations or not
       miinn_nsa,&                  ! Max # of iterations
       mcfll_nsa,&                  ! Max # of CFL values for the CFL_LOCAL option
       miinn_pseud_nsa,&            ! Max # pseudo-time of iterations
       last_iters_nsa,&             ! Solver iterations for the current sub-iteration
       msoit_nsa,&                  ! Max # of solver iterations
       npica_nsa,&                  ! Number of Picard iteration (Newton's lin.)
       kranr_nsa(3,4),&             ! Rayleigh ranges type (3 space dimensions, lower and greater)
       nkryd_nsa,&                  ! Krylov dimension
       neule_nsa,&                  ! # of Euler time steps
       nisaf_nsa,&                  ! Initial time step for variable cfl 
       kfl_safet_table_nsa,&        ! Safety factor table flag
       minew_nsa,&                  ! Initial Newton iteration
       nunkn_nsa                    ! # of unknonws = ndofn_nsa*npoin  

  real(rp) ::&
       vranr_nsa(3,3),&             ! Rayleigh ranges limits
       frayl_nsa(3,4),&             ! Rayleigh parameters (3 free parameters, acoustic and gravity)
       dtcri_nsa( 10),&             ! Critical time step
       dtlim_nsa     ,&             ! Factor to limit the local time step
       dtmax_nsa( 10),&             ! Maximum of the local time steps
       shock_nsa    ,&              ! Shock capturing parameter
       shtol_nsa    ,&              ! Shock capturing tolerance (a la alvaro)
       sstol_nsa    ,&              ! Steady state tolerance
       cotol_nsa    ,&              ! Convergence tolerance
       corat_nsa    ,&              ! Convergence tolerance
       penal_nsa    ,&              ! Penalization factor
       dtext_nsa    ,&              ! Externally fixed time step 
       safet_nsa    ,&              ! Safety factor for time step
       safet_table_nsa(2,10) ,&     ! Safety factor table: value, residual
       safex_nsa    ,&              ! Safety factor for time step
       safma_nsa    ,&              ! Safety factor for time step
       safet_pseud_nsa    ,&        ! Safety factor for pseudo time step
       safrk_nsa    ,&              ! Delta time shrinking factor for RK schemes 
       theta_nsa(10),&              ! "Explicitness" weighting factors
       solco_nsa    ,&              ! Solver tolerance
       resid_nsa(4) ,&              ! Residual for outer iterations (Mom,Cont,Ene,Glo)
       resou_nsa(4) ,&              ! Reference residual for outer iterations (Mom,Cont,Ene,Glo)
       resin_nsa(4) ,&              ! Reference residual for inner (sub)iterations (Mom,Cont,Ene,Glo)
       resin_first_nsa(4) ,&        ! Reference residual for first inner (sub)iterations (Mom,Cont,Ene,Glo)
       xfree_nsa,    &              ! X Coordinate of the plane where to free
       weigh_nsa,    &              ! Weight of dU/dt in the residual
       err01_nsa(2),&               ! L1 error u
       err02_nsa(2),&               ! L2 error u
       err0i_nsa(2),&               ! Linf error u
       err11_nsa(2),&               ! L1 error grad(u)
       err12_nsa(2),&               ! L2 error grad(u)
       err1i_nsa(2),&               ! Linf error grad(u)
       fcfll_nsa(500,2),&           ! Local CFL values
       cpu_nstal(2)                 ! CPU for the NST problem

  real(rp), pointer ::&
       unkna_nsa(:),&               ! RK: linear combination of previous iterations unknowns
       unkit_nsa(:,:)               ! RK: previous iterations unknowns
       
!------------------------------------------------------------------------
! Boundary conditions: read, modified or defined  in nsa_reabcs
!------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_confi_nsa,&              ! Confined flow
       kfl_local_nsa,&              ! Local system of reference
       kfl_inlet_nsa,&              ! b.c. to force the mass flow rate at the inlet
       kfl_conbc_nsa,&              ! Constant b.c.
       kfl_dumbo_nsa,&              ! Dump derived boundary conditions
       kfl_spong_nsa,&              ! Rayleigh sponge (old defined by Mariano) Obsolete
       kfl_sponge_nsa,&             ! Sponge flag (SM)
       kfl_skewb_nsa,&              ! Skew symetric (skcos) label in bvess
       kfl_bofty_nsa,&              ! Boundary conditions file type (nstinc or nastal)
       kfl_bofie_nsa,&              ! Boundary conditions fixed from fields, when present
       kfl_nopen_nsa,         &      ! No penetration condition
       kfl_tredg_nsa,         &      ! Free boundary nodes at edges
       nbhie_nsa,&                  ! Boundary condtions hierarchy
       nodpr_nsa,&                  ! Node on which pressure is prescribed
       mtloa_nsa(10), &             ! Number of data points for the transient boundary functions
       nfunc_nsa                    ! Number of transient boundary conditions function
  real(rp) ::&
       delta_nsa,&                  ! Distance to the wall
       inlet_nsa,&                  ! Target mass flow rate at the inlet
       rtico_nsa(20,2), &           ! Time counters for the transient boundary conditions
       angle_tredg_nsa              ! Angle threshold for the trailing edge detection (DEGREES)

  type(bc_nodes), pointer              :: &
       tncod_nsa(:)                          ! Node code type
  type(bc_nodes), pointer              :: &
       tgcod_nsa(:)                          ! Geometrical node code type
  type(bc_bound), pointer              :: &
       tbcod_nsa(:)                          ! Boundary code type
  type(r2p),     pointer                :: &
       tload_nsa(:)                          ! time dependent boundary functions (to be read in a file)


  character(5)   ::                 chhie_nsa(10)

!  real(rp),     allocatable ::&
!       thdif_nsa(:)                 ! Thermal diffusion
!
! Boundary conditions and geometry vectors
! 
  integer(ip),  pointer      ::&
       kfl_fixno_nsa(:,:),&         ! Nodal fixity 
       kfl_fixbo_nsa(:),&           ! Element boundary fixity
       kfl_fixrs_nsa(:),&           ! Reference system for the BV
       kfl_funbo_nsa(:),&           ! Functions for boundary element bc
       kfl_funno_nsa(:)             ! Functions for node bc
  real(rp),     pointer      ::&
       rhsou_nsa(:)    ,&           ! Residual for the Newton iteration that works as a source        
       reafo_nsa(:)    ,&           ! Residual Ax - b
       jacrot_du_dq_nsa(:,:,:),&                     ! Jacobian+Rotation matrix  Q -> U  (ANTES ROTQU)
       jacrot_dq_du_nsa(:,:,:),&           ! Transpose of the Jacobian+Rotation matrix  U -> Q (ANTES ROTUQ)
       bvnat_nsa(:,:),&             ! Natural bc values
       bvess_nsa(:,:,:),&           ! Essential bc (or initial) values
       rekee_nsa(:,:),  &           ! Hydrostatic reference field. Order: velocity components, density, temperature, pressure and Brunt
       densiref_nsa(:),  &         !  Reference field read from sounding
       temperef_nsa(:),  &         !  Reference field read from sounding
       qvaporef_nsa(:),  &         !  Reference field read from sounding
       uveloref_nsa(:),  &         !  Reference field read from sounding
       vveloref_nsa(:),  &         !  Reference field read from sounding
       bspon_nsa(:,:),  &          ! Rayleigh sponge coefficients vector
       skcos_nsa(:,:,:)            ! Cosine matrices of NS skew systems

  !Column-wise variables for kessler
  real(rp),     pointer  ::&
       rhocol_nsa(:),    &!
       tcol_nsa(:),      &!
       qvcol_nsa(:),     &!
       qccol_nsa(:),     &!
       qrcol_nsa(:),     &!
       vtcol_nsa(:),     &!
       prodcol_nsa(:),   &!
       prodkcol_nsa(:),  &!
       vtdencol_nsa(:),  &!
       rdzkcol_nsa(:),   &!
       rdzwcol_nsa(:),   &!
       rhokcol_nsa(:),   &!
       factorcol_nsa(:), &!
       pcol_nsa(:),      &!
       z_nsa(:),         &!
       z_col_nsa(:),     &!
       q_nsa(:,:),       &!
       qref_nsa(:,:),    &!
       q_col_nsa(:,:),   &!
       qref_col_nsa(:,:),&!
       rainnc_nsa(:),    &!
       rainncv_nsa(:),   &!
       prodk_nsa(:),     &!
       aa_nsa(:),        & !aa and bb: Sponge coefficients
       bb_nsa(:)

  integer(ip)                ::&
       kfl_funty_nsa(20,20)           ! Function type and number of paremeters            
  real(rp)                   ::&
       fubcs_nsa(20,20),       &
       zspon_nsa(3)                   ! Sponge heights    
  type(r1p),    pointer      ::&
       funpa_nsa(:)                 ! Function parameters

!------------------------------------------------------------------------
! Output and Postprocess: read, modified or defined  in nsa_reaous
!------------------------------------------------------------------------
  integer(ip)                ::&
       kfl_bodyf_nsa,&              ! Body aerodynamic coefficients flag
       kfl_pro2d_nsa,&              ! 2D profile distribution 
       kfl_exacs_nsa,&              ! Exact solution for the NS eq.
       kfl_crens_nsa,&              ! Compute or not vector crens (master local)
       kfl_chkpo_nsa(2)             ! Checkpoint-restart file flags: input(1) and output(2)
 
  real(rp) ::&
       bodyr_nsa(10,2,3),&          ! Body aerodynamics coefficients bounding boxes
       avtim_nsa                    ! Accumulated time for time-averaging
 
  character(150)         :: &
       fil_chkpo_nsa(2)             ! Checkpoint-restart file names: input(1) and output(2)
!--END REA GROUP
!------------------------------------------------------------------------
! Module fields
!------------------------------------------------------------------------
  real(rp),     pointer      ::&
       umoss_nsa(:,:,:),denss_nsa(:,:),eness_nsa(:,:),&        ! Nodal subscale field for U, rho, E
       umosg_nsa(:,:,:,:),densg_nsa(:,:,:),enesg_nsa(:,:,:),&    ! Gauss subscale field for U, rho, E
       ortpr_nsa(:,:),&       ! Nodal orthogonal projection onto the FE space W_h
       shecp_nsa(:),&         ! Specificic heat Cp 
       avpre_nsa(:),&         ! Time-averaged pressure
       avtem_nsa(:),&         ! Time-averaged temperature
       avvel_nsa(:,:),&       ! Time-averaged velocity
       avve2_nsa(:,:),&       ! Time-averaged normal stress 
       avmom_nsa(:,:),&       ! Time-averaged momentum 
       avvxy_nsa(:,:)         ! Time-averaged shear stress

  type(r3p), pointer ::&
       shocktau_nsa(:)        ! Shock capturing diffusion per element / gauss


!------------------------------------------------------------------------
! Derived variables (slave-local in MPI)
!------------------------------------------------------------------------
  logical                  :: &
       iloner,weparal
  integer(ip)              :: &
       izone_nsa ,&                 ! Zone for Nastal
       nbdfp_nsa ,&                 ! Number of terms in temporal derivative
       ndtdf_nsa ,&                 ! Time increment dof
       nindx_nsa(3,3) ,&            ! Index for the hessian (voigt -> cartesian)
       kfl_tiaor_nsa,&              ! Original time accuracy
       kfl_refre_nsa(6) ,&          ! Compute the reference residual for the nsa_cvgunk's itask  
       ipose_nsa                    ! Postprocess counter (Ensight)
  real(rp)                 :: &
       stapa_nsa(10) ,&             ! Auxiliar vector for state law parameters
       pabdf_nsa(10)                ! BDF parameters
  real(rp)                 :: &
       parkb_nsa(8,2)    ,&         ! Time integration parameters for 
       parka_nsa(8,8,2)           !   high order time schemes: A matrix and b vector

  type(i1p),   pointer     :: &
       leobl_nsa(:)                 ! Anti lboel: el-2-bo correspondence 

  real(rp),    pointer     ::&
       cosma_nsa(:,:,:,:,:,:) ,&      ! Fourier space stabilization coefficients (cosine)
       sinma_nsa(:,:,:,:,:,:) ,&      ! Fourier space stabilization coefficients (sine)
       dtieq_nsa(:,:,:),&           ! Dt per node, per equation
       tauti_nsa(:,:),&             ! Taus
       vmacp_nsa(:),&               ! Lumped mass matrix for a macro node
       vdiag_nsa(:),&               ! Diagonal implicit matrix
       dunkn_nsa(:),&               ! Delta unknown
       resmo_nsa(:),&               ! Smoothed residual
       rhsax_nsa(:,:,:)             ! Temporal internal residuals
  integer(ip),    pointer     ::&
       linod_nsa(:),liele_nsa(:)    ! Auxiliar node or element lists    

!------------------------------------------------------------------------
! Global parameters: drag and lift coefficients, forces, pitches, etc.
!------------------------------------------------------------------------

  real(rp)                ::&
       afact_nsa,           &      ! iso_wrf A factor
       vemxm_nsa(2),        &      ! max-min values for the velocity module
       vamxm_nsa(2,10)             ! Max-min values for the variables
  integer(ip)             ::&
       kfl_zevel_nsa(2)            ! zero velocity at the start

!
! Aerodynamic body forces 
!  
  real(rp),    target     :: &
       ccoef_nsa(20)                 ! coefficients and forces

  real(rp)                ::&
       fbody_nsa(2,3),surto_nsa,sufro_nsa(3),&
       clise_nsa(msets_nsa),cdrse_nsa(msets_nsa),&
       sforc_nsa(msets_nsa,3),sptch_nsa(msets_nsa,3)
  real(rp),    pointer     ::&
       clicu_nsa(:),cdrcu_nsa(:)              ! Lift and drag per cut
  integer(ip)             ::&
       npopr_nsa

!------------------------------------------------------------------------
! Convection-Diffusion-Reaction Equation
!------------------------------------------------------------------------
  real(rp)                :: diffu_nsa ,react_nsa
  real(rp)                :: conve_nsa(3)

    !>-----------------------------------------------------------------------<!
    ! LODI
    !>-----------------------------------------------------------------------<! 
    logical(ip)       :: euler_nsa 
!    real(rp), pointer ::  nsa_forfsi(:,:)       !> force for fsi 
    real(rp), pointer :: xschlrn_nsa(:,:)      !> schlieren 
   !real(rp), pointer ::  schlrn_nsa(:)        !> schlieren 
   !real(rp), pointer ::   xsrc_nsa(:,:,:)     !> gauss point source 
    real(rp), pointer ::  xchrc_nsa(:,:,:,:)   !> gauss point characteristic 
   !real(rp), pointer ::   chrc_nsa(:,:,:)     !> nodal point characteristic 
   !real(rp), pointer ::  sxinv_nsa(:,:,:,:,:) !> matrix S**-1
    real(rp), pointer ::     sx_nsa(:,:,:,:,:) !> matrix S
    real(rp), pointer ::  vmach_nsa(:)         !> nodal point Mach number 
    real(rp), pointer ::  gamma_nsa(:)         !> nodal point gamma  
    real(rp), pointer ::  sound_nsa(:)         !> nodal point gamma    
    real(rp), pointer ::  xdvol_nsa(:,:)       !> (nelem,mgaus)
    real(rp)          :: sigma_lodi_nsa        !> sigma_lodi_nsa = sigma/Lx 
    real(rp)          :: kfact_lodi_nsa        !> kfact = sigma_lodi_nsa * (1-Mach**2)
    real(rp)          :: prefe_lodi_nsa        !> Kfact(Pi - Pref) 
    integer(ip)       :: lodi_nsa              !> actived 

    !> nsa_lodi_normaltype   
    type lodi_type
      integer(ip) :: id, idime, ichrc
    endtype lodi_type

    type(lodi_type), pointer :: normal_nsa(:)
    !>-----------------------------------------------------------------------<! 
contains
  
  function nsa_zonenode(ipoin)
    use def_kintyp
    use def_domain
    use def_master
    implicit none
    
    logical :: nsa_zonenode
    integer(ip) :: izone,ipoin
    
    nsa_zonenode = .true.
    izone= lzone(ID_NASTAL)
    
    if (nzone > 1) then
       !    if (lpoiz(izone,ipoin) == 0) nsa_zonenode = .false.
    end if
    
  end function nsa_zonenode
    

end module def_nastal
