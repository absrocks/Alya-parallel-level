module def_radiat
  !------------------------------------------------------------------------
  !****f* Radiat/def_radiat
  ! NAME 
  !    def_radiat
  ! DESCRIPTION
  !    Heading for the Radiat routines
  ! USES
  ! USED BY
  !    Almost all
  !***
  !------------------------------------------------------------------------
  use def_kintyp

  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip), parameter ::&
       lun_bound_rad = 210, lun_psmat_rad = 214, lun_funck_rad = 215, &
       lun_funcc_rad = 216, lun_intbc_rad = 221, lun_dynin_rad = 222, &
       lun_dynou_rad = 223, lun_dynlo_rad = 224, lun_dynre_rad = 225, & 
       lun_splot_rad = 232, lun_ramsh_rad = 233, lun_rares_rad = 234
       
  character(150)                        :: &
       fil_ramsh_rad,                      &
       fil_rares_rad,                      &
       fil_dynin_rad,                      &  ! Input dynamic model
       fil_dynou_rad                          ! Output dynamic model

  real(rp),      parameter :: &
       zetem = epsilon(1.0_rp),            &
       Steph_rad =  0.00000005670373_rp            ! Stephan-Boltzmann constant

  integer(ip),   parameter              :: &
       nexap_rad=10                         ! # exact solution parameters

  !------------------------------------------------------------------------
  ! Physical problem: read in rad_reaphy
  !------------------------------------------------------------------------

  !
  ! Physical properties / Materials
  !
!--BEGIN REA GROUP
  real(rp),     allocatable             :: &
       scatt_rad(:),                   & ! Scattering coeff sigma_S
       aniso_rad(:),                   & ! Anisotropy coeff C
       absor_rad(:)                      ! Absorption coeff a  
    
  integer(ip)                           :: &
       nspec_rad,                      & ! Number of species in play
       kfl_parti_rad,                 & ! Are there particles in suspension 
       idtem_rad                         ! ID of particles temperature
 
!!$  type(adrtyp), target                  :: &
!!$       adreq_rad(1)                          ! ADR eqn: type  !!F Do I need this?

  !------------------------------------------------------------------------
  ! Numerical problem: read in rad_reanut
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_ellen_rad,                      & ! =0,1 for min/max element length
       kfl_sgsno_rad,                      & ! Subscale non-linear tracking
       kfl_taust_rad,                      & ! Tau calculation option
       kfl_ortho_rad,                      & ! Orthogonal SGS
       kfl_limit_rad,                      & ! Limiter
       kfl_bubbl_rad,                      & ! Bubble
       kfl_assem_rad,                      & ! Assembly
       kfl_normc_rad,                      & ! Norm of convergence
       miinn_rad,                          & ! Max # of iterations
       misgs_rad                             ! Max # of SGS iterations

  real(rp)                              :: &
       staco_rad(3),                       & ! Stability constants
       cotol_rad,                          & ! Convergence tolerance
       relsg_rad,                          & ! Relaxation parameter of subgrid scale
       tosgs_rad                             ! Subgrid scale tolerance

  !------------------------------------------------------------------------
  ! Output and Postprocess: read in rad_reaous
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_splot_rad,                      & ! Output for 3D gnuplot flag
       kfl_psmat_rad,                      & ! PS file of matrix profile
       kfl_exacs_rad,                      & ! Exact solution for the heat eq.
       kfl_atest_rad,                      & ! Testing functions for temperature and density
       npp_bound_rad                         ! Postprocess boundary conditions
  real(rp)                              :: &
       expar_rad(nexap_rad)                  ! Exact or test solution parameters


  !------------------------------------------------------------------------
  ! Boundary conditions: read in rad_reabcs
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_conbc_rad,                      & ! Constant boundary conditions
       kfl_intbc_rad,                      & ! Interpolate boundary conditions
       kfl_inico_rad,                      & ! Initial conditions
       npnat_rad                             ! # variable natural bc
  real(rp)                              :: &
       bvcoe_rad(10)                         ! Initial condition coefficients
  integer(ip),   pointer                :: &
       kfl_funty_rad(:)                      ! Function type  
  real(rp),      pointer                :: &
       funpa_rad(:,:)                        ! Function parameters
  type(bc_nodes), pointer               :: &     
       tncod_rad(:)                          ! Node code type
  type(bc_nodes), pointer               :: &     
       tgcod_rad(:)                          ! Geometrical node code type
  type(bc_bound), pointer               :: &     
       tbcod_rad(:)                          ! Boundary code type
!--END REA GROUP
  !------------------------------------------------------------------------
  ! Boundary conditions
  !------------------------------------------------------------------------

  integer(ip),   pointer                :: &
       kfl_fixno_rad(:,:),                 & ! Nodal fixity 
       kfl_fixbo_rad(:),                   & ! Element boundary fixity
       kfl_funno_rad(:),                   & ! Functions for node bc       
       kfl_funbo_rad(:)                      ! Functions for boundary bc       
  real(rp),      pointer                :: &
       bvess_rad(:,:),                     & ! Essential bc values
       bvnat_rad(:,:,:)                      ! Natural bc values

  !------------------------------------------------------------------------
  ! Bubble treatment
  !------------------------------------------------------------------------

  real(rp) ::&  
       bumat_rad,                          & ! Bubble matrix coefficient
       burhs_rad                             ! Bubble RHS
  real(rp),     allocatable             :: &
       rtemp_rad(:,:),                     & ! Residual at parent element Gauss points
       ttemp_rad(:,:),                     & ! Adjoint at parent element Gauss points
       rtemb_rad(:,:,:),                   & ! Residual at mini-element Gauss points
       ttemb_rad(:,:,:),                   & ! Adjoint at mini-element Gauss points
       rcoef_rad(:,:)                        ! Residual coefficients
  real(rp),     pointer                 :: & 
       fixnb_rad(:)                          ! Bubble Fixity
  integer(ip),  pointer                 :: &
       lmatb_rad(:)                          ! Bubble material

  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       ncomp_rad,                          & ! Number of components of the radiation
       kfl_sgsve_rad,                      & ! If velocity subgrid scale is considered
       kfl_goite_rad,                      & ! Keep iterating
       nzmat_rad,                          & ! Matrix size
       nzrhs_rad,                          & ! RHS size
       iknbo_rad,                          & ! =1 if boundary number not known
       ittot_rad                             ! Total number of iterations
  real(rp)                              :: &
       dtinv_rad,                          & ! 1/dt
       resid_rad,                          & ! Residual for outer iterations
       pabdf_rad(10),                      & ! BDF factors
       ramin_rad,                          & ! Minimum radiation intensity
       ramax_rad,                          & ! Maximum radiation intensity
       err01_rad(2),                       & ! L1 error T
       err02_rad(2),                       & ! L2 error T
       err0i_rad(2),                       & ! Linf error T
       err11_rad(2),                       & ! L1 error grad(T)
       err12_rad(2),                       & ! L2 error grad(T)
       err1i_rad(2)                          ! Linf error grad(T)
  real(rp), target                      :: &
       resgs_rad(2)                          ! SGS residual
  type(r1p),    allocatable             :: &
       viewf_rad(:)                          ! View factors Fij (radiation)
  real(rp),     pointer                 :: &
       tepro_rad(:)                          ! Radiation orthogonal projection
  ! 
  ! Variables for radiation
  !
  real(rp),     pointer                 :: &
       radav_rad(:,:),                     & ! Average radiation intensity
       grrad_rad(:,:),                     & ! Radiation gradients
       raold_rad(:),                       & ! Old radiation
       conce_rad(:,:),                     & ! Concentration of species for use in radiation
       tempe_rad(:)                          ! Temperature for use in radiation

  type(r2p),   pointer     :: &
       rasgs_rad(:)                         ! Radiation subgrid scale    RADIAT

end module def_radiat
