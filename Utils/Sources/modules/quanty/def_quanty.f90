module def_quanty
  !------------------------------------------------------------------------
  !****f* Quanty/def_quanty
  ! NAME 
  !    def_quanty
  ! DESCRIPTION
  !    Heading for the Quanty routines
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
       lun_pdata_qua = 1501, lun_outpu_qua = 1502, lun_conve_qua = 1503,&
       lun_solve_qua = 1504, lun_rstar_qua = 1505,                      &
       lun_witne_qua = 1506, lun_ppseu_qua = 1507       

  character(150)                        :: &
       fil_rstar_qua,                      &
       fil_conve_qua,                      &
       fil_pdata_qua,                      &
       fil_solve_qua,                      &
       fil_outpu_qua,              &
       fil_witne_qua,                      &
       fil_ppseu_qua

  real(rp),      parameter :: &
       zequa = epsilon(1.0_rp)
  integer(ip),   parameter              :: &
       nvars_qua=10,                       &  ! # set variables
       nvarp_qua=20,                       &  ! # postprocess variables
       nvart_qua=10,                       &  ! # postprocess times 
       nvarw_qua=10,                       &  ! # witness point
       ncoef_qua=10,                       &  ! # coefficient for properties
       nsetp_qua=10,                       &  ! # set parameters
       mwitn_qua=100                          ! Max # witness points

  !------------------------------------------------------------------------
  ! Physical problem: read in qua_reaphy
  !------------------------------------------------------------------------
!--BEGIN REA GROUP
  integer(ip)                           ::  &
       kfl_timei_qua,                       & ! temporal problem d()/dt
       kfl_dftgs_qua,                       & ! DFT ground state
       kfl_potxc_qua,                       & ! XC potential type
       kfl_alele_qua,                       & ! all electron problem
       kfl_nolin_qua,                       & ! nonlinar problem
       kfl_perio_qua,                       & ! periodic problem
       kfl_coulo_qua,                       & ! Coulomb potential
       kfl_bfieldx_qua,                     &
       kfl_bfieldy_qua,                     &
       kfl_bfieldz_qua,                     &
       klf_btemp_law,                       & ! B temp. law
       kfl_efieldx_qua,                     &
       kfl_efieldy_qua,                     &
       kfl_efieldz_qua,                     &
       klf_etemp_law,                       & ! E temp. law
       kfl_spinb_qua,                       & 
       kfl_spiorb_qua,                      &
       kfl_relat_qua,                       &
       kfl_vother_qua,                      &
       ncuanpal_qua,                        &
       lcuanorb_qua,                        &
       nspin_qua,                           &
       lawma_qua,                           &
       law_vother_qua,                      &                              
       natoms_qua,                          & !numero de atomos
       nespecies_qua,                       & !numero de especies
       ncomp_eig,                           & ! Number of autoestates
       nestates,                            & ! Number of estates to be calculed in DFT case
       noutput                                ! shicht de salida

  real(rp)                               :: &
       eig_evol_qua,                        &  !autovalor que evoluciona
       coulo_qua,                           &  ! Value of Coulomb potential
       bfieldx_qua,                         &
       bfieldy_qua,                         &
       bfieldz_qua,                         &
       btemp_law_qua,                       & ! B temp. law
       efieldx_qua,                         &
       efieldy_qua,                         &
       efieldz_qua,                         &
       etemp_law_qua,                       & ! E temp. law
       frecuencie,                          & ! frecuencia de la ley temporal 
       spinb_qua,                           &
       spiorb_qua,                          &
       relat_qua,                           &
       vother_qua,                          &
       w_vother_qua,                        &
       x0_qua,                              &
       y0_qua,                              &
       z0_qua,                              &
       massa_qua,                           &
       mezcla                                 ! mezcla de densidades en DFT 

  type(atomo), pointer                   :: &
       atomo_qua(:) 

  type(especie), pointer                 :: &
       especie_qua(:) 

  integer(ip), allocatable               ::  noccupa(:)   ! Ocupation by nestates

  !------------------------------------------------------------------------
  ! Numerical problem: read in qua_reanut
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       miinn_qua,                          & ! Max # of iterations
       kfl_tiacc_qua,                      & ! Time accuracy
       kfl_tisch_qua                         ! ime integration scheme
  real(rp)                              :: &
       cotol_qua                             ! Convergence tolerance

  type(soltyp), pointer                 :: &
       solve_qua(:)                          ! Solver type
  type(eigtyp), pointer                 :: &
       eigen_qua(:)                         ! Eigen Solver type

  !------------------------------------------------------------------------
  ! Output and Postprocess: read in qua_reaous
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       npp_inits_qua,                      & ! Postprocess initial step
       npp_iniso_qua,                      & ! Postprocess initial solution
       npp_stepi_qua(nvarp_qua),           & ! Postprocess step interval
       npp_setse_qua(nvars_qua),           & ! Postprocess element sets calculation
       npp_setsb_qua(nvars_qua),           & ! Postprocess boundary sets calculation
       npp_setsn_qua(nvars_qua),           & ! Postprocess node sets calculation
       npp_witne_qua(nvarw_qua),           & ! ?
       kfl_exacs_qua                         ! Exact solution for the heat eq.

  real(rp)                              :: &
       pos_tinit_qua,                      & ! Postprocess initial time
       pos_times_qua(nvart_qua,nvarp_qua), & ! Postprocess times       
       setpa_qua(nsetp_qua)                  ! Sets parameters

  !------------------------------------------------------------------------
  ! Boundary conditions: read in qua_reabcs
  !------------------------------------------------------------------------

  integer(ip),   pointer                :: &
       kfl_fixno_qua(:,:),                 & ! Nodal fixity 
       kfl_fixbo_qua(:),                   & ! Element boundary fixity
       kfl_funno_qua(:),                   & ! Functions for node bc       
       kfl_funbo_qua(:),                   & ! Functions for boundary bc       
       kfl_funty_qua(:)                      ! Function type 

  real(rp),      pointer                :: &
       bvess_qua(:,:),                     & ! Essential bc values
       bvessH_qua(:,:),                    & ! Essential bc values for Poisson
       funpa_qua(:,:)                        ! Function parameters

  integer(ip)                           :: &
       kfl_conbc_qua,                      & ! Constant boundary conditions
       kfl_inidi_qua,                      & ! Initial problem
       kfl_inico_qua,                      & ! Initial conditions
       kfl_intbc_qua,                      & ! Initial conditions
       npnat_qua                             ! # variable natural bc
!--BEGIN REA GROUP
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  character(5)                          :: & 
       wopos_qua(2,nvarp_qua)                ! Name and character of the postprocess variables
  character(5)                          :: & 
       woese_qua(nvars_qua),               & ! Name and character of the element set variables
       wobse_qua(nvars_qua),               & ! Name and character of the boundary set variables
       wonse_qua(nvars_qua),               & ! Name and character of the node set variables
       wowit_qua(nvarw_qua)                  ! Name and character of the witness point variables
  type(r1p),    allocatable             :: &
       viewf_qua(:)                          ! View factors Fij (radiation)
  real(rp),     pointer                 :: &
       veset_qua(:,:),                     & ! Set element values
       vbset_qua(:,:),                     & ! Set boundary values
       vnset_qua(:,:)                        ! Set node values

  ! 
  ! Others
  !
  real(rp),     pointer                 :: &
       colec(:)                              ! colec !

  real(rp), allocatable               ::  v_pot_ps(:)   ! potencial construido con PS 
  real(rp), allocatable               ::  v_xc(:)       ! potencial construido con XC 
  real(rp), allocatable               ::  v_hartree(:)  ! potencial construido con hartree 
  complex(rp) , ALLOCATABLE           ::  NonLoc(:,:,:)
  real(rp), allocatable               ::  NL_denom(:,:)

  integer(ip)                           :: &
       ncomp_qua,                         & ! Number iter temp
       kfl_stead_qua,                     & ! Steady-state has been reached
       kfl_tiaor_qua,                     & ! Original time accuracy
       kfl_goite_qua,                     & ! keep iteratin
       ittot_qua,                         & 
       lmaximo                              ! maximo l in PPs

  real(rp)                              :: &
       sstol_qua,                          & !tol for time 
       dtinv_qua,                          & ! 1/dt
       dtcri_qua,                          & ! Critical time step
       resid_qua,                          & ! Residual for outer iterations
       pabdf_qua(10),                      & ! BDF factors
       rhomin_qua,                          & ! Minimum temperature
       rhomax_qua,                          & ! Maximum temperature
       err01_qua(2),                       & ! L1 error T
       err02_qua(2),                       & ! L2 error T
       err0i_qua(2),                       & ! Linf error T
       err11_qua(2),                       & ! L1 error grad(T)
       err12_qua(2),                       & ! L2 error grad(T)
       err1i_qua(2),                       & ! Linf error grad(T)
       avtim_qua                             ! Averaging time

end module def_quanty
