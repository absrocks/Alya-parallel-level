!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_reanut.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Porous specific variables
!> @details Porous specific variables
!> @} 
!------------------------------------------------------------------------
module def_porous
  use def_kintyp

  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip), parameter ::&
       lun_bound_por = 2410, lun_psmat_por = 2414, lun_funck_por = 2415, &
       lun_funcc_por = 2416, lun_intbc_por = 2421, lun_dynin_por = 2422, &
       lun_dynou_por = 2423, lun_dynlo_por = 2424, lun_dynre_por = 2425, & 
       lun_splot_por = 2432, lun_ramsh_por = 2433, lun_rares_por = 2434

  !  character(150)                        :: &
  !       fil_ramsh_por,                      &
  !       fil_rares_por

  real(rp),      parameter :: &
       zepor = epsilon(1.0_rp)
  integer(ip),   parameter              :: &
       nprsa_por=2,                        & ! 2 options pressure or saturation
       ndels_por=1,                        & ! Number of steps to delay saturation solution
       mwell_por=100,                      & ! maximum number of wells + injectors
       kfl_malta_por=1                       ! Option from paper from malta being tested ! OJO /=0 solo listo para g=0
  ! ademas solo esta listo para el caso en que damos pbh o Qw - para los casos Qt aun 
  ! queda por analizar escribir y luego programar - Tambien falta pbh en inyector.
  !
  ! Type
  !
  type well_type
     integer(ip)          :: itype  
     integer(ip)          :: ntime 
     integer(ip)          :: nposi 
     real(rp)             :: radiu 
     real(rp),    pointer :: q_table(:,:)
     real(rp),    pointer :: pbh_table(:,:)
     real(rp),    pointer :: pbh_coord(:)     
  end type well_type
  type(well_type), pointer :: tywel_por(:)

  !------------------------------------------------------------------------
  ! Physical problem: read in por_reaphy
  !------------------------------------------------------------------------
  !--BEGIN REA GROUP
  integer(ip)                           :: &
       kfl_timei_por(nprsa_por),           & ! Existence of dT/dt
       kfl_advec_por(nprsa_por),           & ! Existence of convective term in pressure or saturation equation
       nmate_por,                          & ! Number of materials for tkrel_por
       nwell_por,                          & ! Number of wells
       mrows_por,                          & ! Maximum number of rows for tkrel_por
       mroww_por,                          & ! Maximum number of rows for wells botom hole pressure
       mheiw_por                             ! maximum number of points per well

  integer(ip),  pointer                 :: &  
       nrows_por(:),                       & ! number of rows for each material in tkrel_por
       nroww_por(:),                       & ! number of rows for the Pbh of each well 
       nheiw_por(:),                       & ! number vertical points that form the well -- read from data file 
       kfl_wellc_por(:)                      ! Well condition ( 1 Prescribed water flow rate , 2 Presc. Total  ....)

  integer(ip),  pointer                 :: & ! Using fields - sendat not needed
       iwell_por(:)                          ! number of well to which ipoin belongs

  real(rp)                  :: &
       comro_por    ,          &      ! Rock compressibility
       comwa_por    ,          &      ! Water compressibility
       comoi_por    ,          &      ! Oil compressibility
       bwref_por    ,          &      ! Bw at reference Pressure
       boref_por    ,          &      ! Bo at reference Pressure
       prref_por    ,          &      ! Reference Pressure
       muwat_por    ,          &      ! Water viscosity
       muoil_por    ,          &      ! Oil viscosity
       denwa_por    ,          &      ! Water density at reference pressure
       denoi_por    ,          &      ! Oil density at reference pressure
       denhy_por    ,          &      ! Density used to substract Hydrostatic component
       gravi_por(3) ,          &      ! Gravity vector
       grnor_por    ,          &      ! Gravity norm
       prini_por                      ! Initial constant Pressure

  real(rp),    pointer      :: &
       tkrel_por(:,:,:),       &        ! Table for k_rw & k_ro
       wvalu_por(:,:,:),       &        ! Well values: it can be the bottom hole pressure or flow rate table
       rwell_por(:)                     ! Well radius

  real(rp),    pointer      :: &        ! Using fields - sendat not needed
       poro0_por(:),           &        ! Porosity at reference pressure
       perme_por(:,:),         &        ! Permeability
       satin_por(:)                     ! Elemental Initial water saturation

  !------------------------------------------------------------------------
  ! Numerical problem: read in por_reanut
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_ellen_por,                      & ! =0,1 for min/max element length
       kfl_sgsti_por(nprsa_por),           & ! Subscale time tracking
       kfl_sgsno_por(nprsa_por),           & ! Subscale non-linear tracking
       kfl_taust_por(nprsa_por),           & ! Tau calculation option
       kfl_ortho_por(nprsa_por),           & ! Orthogonal SGS
       kfl_limit_por(nprsa_por),           & ! Limiter
       kfl_shock_por(nprsa_por),           & ! Shock capturing type 
       kfl_difwe_por,                      & ! Diffusion to be added at the wells
       kfl_tiacc_por(nprsa_por),           & ! Temporal accuracy
       neule_por,                          & ! Number of Euler time steps 
       kfl_tisch_por,                      & ! Time integration scheme
       kfl_normc_por,                      & ! Norm of convergence
       miinn_por                             ! Max # of iterations
  !       misgs_por                             ! Max # of SGS iterations  ! not used for the moment

  real(rp)                              :: &
       staco_por(3),                       & ! Stability constants
       shock_por,                          & ! Shock capturing parameter
       difwe_por,                          & ! Diffusion to be added at the wells
       safet_por,                          & ! Safety factor for time step
       sstol_por,                          & ! Steady state tolerance
       cotol_por,                          & ! Convergence tolerance
       relax_por,                          & ! Relaxation factor
       bemol_por,                          & ! Integration by parts of convective term
       relsg_por,                          & ! Relaxation parameter of subgrid scale
       tosgs_por                             ! Subgrid scale tolerance

  !------------------------------------------------------------------------
  ! Output and Postprocess: read in por_reabcs
  !------------------------------------------------------------------------

  type(bc_nodes), pointer               :: &     
       tncod_por(:)                          ! Node code type

  !------------------------------------------------------------------------
  ! Output and Postprocess: read in por_reaous
  !------------------------------------------------------------------------

  !  integer(ip)                           :: &
  !       kfl_splot_por,                      & ! Output for 3D gnuplot flag
  !       kfl_psmat_por,                      & ! PS file of matrix profile
  !       kfl_exacs_por,                      & ! Exact solution for the heat eq.
  !       npp_bound_por                         ! Postprocess boundary conditions

  !--END REA GROUP
  !------------------------------------------------------------------------
  ! Boundary conditions
  !------------------------------------------------------------------------

  integer(ip),   pointer                :: &
       kfl_fixno_por(:,:,:)                  ! Nodal fixity 
  real(rp),      pointer                :: &
       bvess_por(:,:,:)                      ! Essential bc values

  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       ncomp_por,                          & ! Number of components of the pressure & saturation
       kfl_sgsve_por,                      & ! If velocity subgrid scale is considered
       kfl_goite_por,                      & ! Keep iterating
       kfl_stead_por,                      & ! Steady-state has been reached 
       nzmat_por,                          & ! Matrix size
       nzrhs_por,                          & ! RHS size
       ittot_por                             ! Total number of iterations
  real(rp)                              :: &
       dtinv_por,                          & ! 1/dt
       dtcri_por,                          & ! Critical time step
       resid_por(2),                       & ! Residual for outer iterations
       pabdf_por(10)                         ! BDF factors


  real(rp), target                      :: &
       resgs_por(2)                          ! SGS residual

  integer(ip),  pointer                 :: &  
       ipwel_por(:,:),                     & ! inverse of iwell_por: for each well and height find ipoin (at each subdomain)
       iheip_por(:)                          ! height level of ipoin


  real(rp),     pointer                 :: &
       nodpo_por(:),                       & ! Nodal projection of the porosity at reference pressure
       nodpe_por(:,:),                     & ! Nodal permebilty, for the moment only for postprocess
       winde_por(:),                       & ! Well index
       wmass_por(:),                       & ! Well mass
       pboho_por(:),                       & ! Well Bottom hole pressure
       dataw_por(:,:),                     & ! Well pressure or flow rate at different heights
       xwell_por(:)

  integer(ip)                           :: &
       kprsa_por                             ! Are we solving pressure(1) or saturation(2) equation

end module def_porous
