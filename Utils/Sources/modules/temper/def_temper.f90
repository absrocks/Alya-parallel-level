module def_temper
  !------------------------------------------------------------------------
  !****f* Temper/def_temper
  ! NAME 
  !    def_temper
  ! DESCRIPTION
  !    Heading for the Temper routines
  ! USES
  ! USED BY
  !    Almost all
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use mod_ADR,   only : ADR_typ

  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip), parameter ::&
       lun_bound_tem = 210, lun_psmat_tem = 214, lun_funck_tem = 215, &
       lun_funcc_tem = 216, lun_intbc_tem = 221, lun_dynin_tem = 222, &
       lun_dynou_tem = 223, lun_dynlo_tem = 224, lun_dynre_tem = 225, & 
       lun_splot_tem = 232, lun_ramsh_tem = 233, lun_rares_tem = 234, &
       lun_lmach_tem = 235
       
  character(150)                        :: &
       fil_ramsh_tem,                      &
       fil_rares_tem,                      &
       fil_dynin_tem,                      &  ! Input dynamic model
       fil_dynou_tem                          ! Output dynamic model

  real(rp),      parameter :: &
       zetem = epsilon(1.0_rp)
  integer(ip),   parameter              :: &
       ncoef_tem=10,                       &  ! # coefficient for properties
       nexap_tem=10,                       &  ! # exact solution parameters
       mkint_tem=100,                      &  ! Max # of knots for k
       mcint_tem=100                          ! Max # of knots for Cp

  !------------------------------------------------------------------------
  ! Physical problem: read in tem_reaphy
  !------------------------------------------------------------------------
!--BEGIN REA GROUP
  integer(ip)                           :: &
       kfl_timei_tem,                      & ! Existence of dT/dt
       kfl_advec_tem,                      & ! Existence of (u.grad)T
       kfl_joule_tem,                      & ! Existence of Joule effect
       kfl_radia_tem,                      & ! Existence of radiation
       kfl_sourc_tem,                      & ! Existence and type of source term
       kfl_cotur_tem,                      & ! Coupling with a turbulence model
       kfl_tfles_tem,                      & ! Thickened flame model activation
       kfl_adiab_tem,                      & ! Calculation of adiabatic mixing
       kfl_condu_tem,                      & ! Conductivity term
       kfl_exint_tem,                      & ! Interpolation of properties
       kfl_inter_tem,                      & ! Interpolation of arrays
       kfl_dynco_tem,                      & ! Dynamical coupling
       kfl_regim_tem,                      & ! Flow regime: incompressible/compressible 
       kfl_prope_tem,                      & ! properties update strategy for CFI model
       kfl_parti_tem,                      & ! There are particles in suspension
       kfl_flux_tem,                       & ! Flag to activate heat flux from fields as BC          
       idtem_tem                             ! ID of particles temperature

  real(rp)                              :: &
       turbu_tem,                          & ! Turbulence parameters
       prtur_tem,                          & ! Turbulent Prandtl number
       cfi_hmax_tem,                       & ! Maximum enthalpy of the mixture for adiabatic calculation
       cfi_hmin_tem,                       & ! Minimum enthalpy of the mixture for adiabatic calculation 
       scond_tem,                          & ! Surface conductivity coefficient
       react_tem                             ! Reaction term
  real(rp),      pointer                :: &
       heat_flux(:,:)                        ! Heat flux prescribed from fields

  !------------------------------------------------------------------------
  ! Numerical problem: read in tem_reanut
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_dttyp_tem,                      & ! Local time step strategy
       kfl_ellen_tem,                      & ! =0,1 for min/max element length
       kfl_sgsti_tem,                      & ! Subscale time tracking
       kfl_sgsno_tem,                      & ! Subscale non-linear tracking
       kfl_taust_tem,                      & ! Tau calculation option
       kfl_ortho_tem,                      & ! Orthogonal SGS
       kfl_limit_tem,                      & ! Limiter
       kfl_shock_tem,                      & ! Shock capturing type 
       kfl_tiacc_tem,                      & ! Temporal accuracy
       kfl_tibub_tem,                      & ! Time integration of bubble
       kfl_plepp_tem,                      & ! PLE activation: =-1 OFF, =0 TEMPER, =3 LOWMACH, =4 ENTHALPY 
       kfl_assem_tem,                      & ! Assembly
       kfl_posit_tem,                      & ! clipping for positive overshoots
       kfl_negat_tem,                      & ! clipping for negative overshoots
       neule_tem,                          & ! Number of Euler time steps 
       kfl_tisch_tem,                      & ! Time integration scheme
       kfl_normc_tem,                      & ! Norm of convergence
       miinn_tem,                          & ! Max # of iterations
       misgs_tem,                          & ! Max # of SGS iterations
       kfl_meshi_tem,                      & ! Mesh interpolator activation flag 
       kfl_discr_tem,                      & ! Discretization method
       kfl_sgsli_tem,                      & ! SGS convection linearization PICARD  
       kfl_sgsac_tem                         ! SGS time accuracy

  real(rp)                              :: &
       staco_tem(3),                       & ! Stability constants
       shock_tem,                          & ! Shock capturing parameter
       safet_tem,                          & ! Safety factor for time step
       source_safet_tem,                   & ! Safety factor for source term contribution to time step
       sstol_tem,                          & ! Steady state tolerance
       cotol_tem,                          & ! Convergence tolerance
       relax_tem,                          & ! Relaxation factor
       bemol_tem,                          & ! Integration by parts of convective term
       relsg_tem,                          & ! Relaxation parameter of subgrid scale
       tosgs_tem,                          & ! Subgrid scale tolerance
       negat_tem,                          & ! negative limit for clipping
       posit_tem                             ! positive limit for clipping

  !------------------------------------------------------------------------
  ! Output and Postprocess: read in tem_reaous
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_splot_tem,                      & ! Output for 3D gnuplot flag
       kfl_psmat_tem,                      & ! PS file of matrix profile
       kfl_exacs_tem,                      & ! Exact solution for the heat eq.
       kfl_viewf_tem(10),                  & ! Output for view factors
       npp_bound_tem                         ! Postprocess boundary conditions
  real(rp)                              :: &
       avtim_tem,                          & ! Averaging Postprocess initial time
       expar_tem(nexap_tem)                  ! Exact solution parameters

       

  !------------------------------------------------------------------------
  ! Boundary conditions: read in tem_reabcs
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_conbc_tem,                      & ! Constant boundary conditions
       kfl_inidi_tem,                      & ! Initial problem
       kfl_inico_tem,                      & ! Initial conditions
       kfl_intbc_tem,                      & ! Initial conditions
       npnat_tem                             ! # variable natural bc
  real(rp)                              :: &
       delta_tem,                          & ! Wall distance
       initial_tem                           ! Initial temperature
  integer(ip),   pointer                :: &
       kfl_funty_tem(:)                      ! Function type  
  real(rp),      pointer                :: &
       funpa_tem(:,:)                        ! Function parameters
  type(bc_nodes), pointer               :: &     
       tncod_tem(:)                          ! Node code type
  type(bc_nodes), pointer               :: &     
       tgcod_tem(:)                          ! Geometrical node code type
  type(bc_bound), pointer               :: &     
       tbcod_tem(:)                          ! Boundary code type
!--END REA GROUP
  !------------------------------------------------------------------------
  ! Boundary conditions
  !------------------------------------------------------------------------

  integer(ip),   pointer                :: &
       kfl_fixno_tem(:,:),                 & ! Nodal fixity 
       kfl_fixbo_tem(:),                   & ! Element boundary fixity
       kfl_funno_tem(:),                   & ! Functions for node bc       
       kfl_funbo_tem(:)                      ! Functions for boundary bc       
  real(rp),      pointer                :: &
       bvess_tem(:,:,:),                   & ! Essential bc values
       bvnat_tem(:,:,:)                      ! Natural bc values

  !------------------------------------------------------------------------
  ! Bubble treatment
  !------------------------------------------------------------------------

  real(rp) ::&  
       bumat_tem,                          & ! Bubble matrix coefficient
       burhs_tem                             ! Bubble RHS
  real(rp),     allocatable             :: &
       rtemp_tem(:,:),                     & ! Residual at parent element Gauss points
       ttemp_tem(:,:),                     & ! Adjoint at parent element Gauss points
       rtemb_tem(:,:,:),                   & ! Residual at mini-element Gauss points
       ttemb_tem(:,:,:),                   & ! Adjoint at mini-element Gauss points
       rcoef_tem(:,:)                        ! Residual coefficients
  real(rp),     pointer                 :: & 
       fixnb_tem(:)                          ! Bubble Fixity
  integer(ip),  pointer                 :: &
       lmatb_tem(:)                          ! Bubble material

  !------------------------------------------------------------------------
  ! Optimization
  !------------------------------------------------------------------------

  real(rp),   pointer ::                   &
       resdiff_tem(:,:),                   & ! Partial Derivative of residual w.r.t design var
       costdiff_tem(:)                       ! Partial Derivative of objective function w.r.t design var
       
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       ncomp_tem,                          & ! Number of components of the temperature
       kfl_sgsve_tem,                      & ! If velocity subgrid scale is considered
       kfl_grdif_tem,                      & ! If there are gradients of conductivity
       kfl_tiaor_tem,                      & ! Original time accuracy
       kfl_goite_tem,                      & ! Keep iterating
       kfl_stead_tem,                      & ! Steady-state has been reached 
       nzmat_tem,                          & ! Matrix size
       nzrhs_tem,                          & ! RHS size
       imate_tem,                          & ! Material under consideration 
       iknbo_tem,                          & ! =1 if boundary number not known
       ittot_tem,                          & ! Total number of iterations
       itsgm_tem,                          & ! maximum number of iterations needed to achieve convergence
       kfl_exist_fixi7_tem,                & ! exists fixity of type 7 
       nunkn_tem                             ! Number of unknowns

  real(rp)                              :: &
       dtinv_tem,                          & ! 1/dt
       dtcri_tem,                          & ! Critical time step
       resid_tem,                          & ! Residual for outer iterations
       pabdf_tem(10),                      & ! BDF factors
       pabds_tem(10),                      & ! BDF factors for sgs
       temin_tem,                          & ! Minimum temperature
       temax_tem,                          & ! Maximum temperature
       err01_tem(2),                       & ! L1 error T
       err02_tem(2),                       & ! L2 error T
       err0i_tem(2),                       & ! Linf error T
       err11_tem(2),                       & ! L1 error grad(T)
       err12_tem(2),                       & ! L2 error grad(T)
       err1i_tem(2),                       & ! Linf error grad(T)      
       vinvt_tem(4),                       & ! inverse of temperature integral, saved each time step
       xmass_tem                             ! total mass, for low Macc
  real(rp), target                      :: &
       resgs_tem(2)                          ! SGS residual
  type(r1p),    allocatable             :: &
       viewf_tem(:)                          ! View factors Fij (radiation)
  real(rp),     pointer                 :: &
       gradc_tem(:,:)                     
  ! 
  ! Others
  !
  real(rp),     pointer                 :: &
       power_tem(:),                       & ! Power of sources (Q)
       avtem_tem(:),                       & ! Average temperature
       avte2_tem(:),                       & ! average tem*tem
       avtev_tem(:,:),                     & ! average vel*tem
       avden_tem(:),                       & ! average density
       avres_tem(:),                       & ! Average residual heat flux
       fvvel_tem(:,:),                     & ! favre average velocity
       grtem_tem(:,:),                     & ! Temperature gradients
       teold_tem(:)                          ! Old temperature

  logical(lg)                           :: &
       kfl_rstar_two                         ! Restarted file was BDF

  integer(ip) :: kfl_code_tem
  !
  ! ADR type
  !
  type(ADR_typ)                         :: &
       ADR_tem                               ! ADR type
 
end module def_temper
