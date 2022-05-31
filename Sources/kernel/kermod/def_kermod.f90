module def_kermod

  !-----------------------------------------------------------------------
  !****f* defmod/def_kermod
  ! NAME
  !   def_kermod
  ! DESCRIPTION
  !   This module is the header of kermod
  !   1 material -> 1 law
  !   1 law      -> 1 responsible module
  !   MATERIAL=1
  !     DENSITY:   LAW=BIFLUID,  PARAMETERS=1.0,2.0,3.0
  !     VISCOSITY: LAW=BIFLUID,  PARAMETERS=1.0,2.0,3.0
  !   END_MATERIAL
  !
  !***
  !-----------------------------------------------------------------------
  
  use def_kintyp_basic,               only : ip,rp,lg,r1p,r2p,r3p
  use def_kintyp_basic,               only : typ_interp
  use def_kintyp_boundary_conditions, only : bc_bound
  use def_kintyp_boundary_conditions, only : bc_nodes
  use def_kintyp_functions,           only : typ_space_time_function
  use def_kintyp_functions,           only : typ_time_function
  use def_kintyp_functions,           only : typ_windk_system
  use def_kintyp_postprocess,         only : witness_geo
  use def_kintyp_postprocess,         only : tyfil
  use def_kintyp_mesh,                only : mesh_type_basic
  use def_coupli,                     only : typ_color_coupling
  use def_domain,                     only : mcodb
  use mod_interp_tab,                 only : typ_tab_coord
  use mod_interp_tab,                 only : typ_lookup_table
  use mod_interp_tab,                 only : typ_lookup_framework
  use mod_interp_tab,                 only : max_lookup_dim

  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  integer(ip),   parameter       ::   &
       mlaws_ker                = 20, & ! Max # laws
       mresp_ker                = 5,  & ! List of responsibles
       mlapa_ker                = 9,  & ! Max # parameters
       mfilt                    = 10, & ! Max # filters
       nelse                    = 20, & ! # Elsest parameters
       max_space_time_function  = 60, & ! Max # space/time functions
       max_time_function        = 20, & ! Max # time functions
       max_windk_systems        = 20, & ! Max # Windkessel systems 
       max_lookup_tab           = 20, & ! Max # lookup tables
       max_lookup_fw            = 10    ! Max # lookup frameworks

  !----------------------------------------------------------------------
  !
  ! Properties
  !
  !----------------------------------------------------------------------
  
  type typ_laws
     integer(ip)                 :: kfl_gradi
     integer(ip)                 :: kfl_deriv
     integer(ip)                 :: kfl_grder
     integer(ip)                 :: kfl_deriv_tur
     integer(ip)                 :: kfl_deriv_vel
     character(5)                :: wname
     integer(ip)                 :: lresp(mresp_ker)
     character(5)                :: where
  end type typ_laws

  type typ_valpr_ker
     integer(ip)                              :: kfl_exist
     integer(ip)                              :: kfl_nedsm
     integer(ip)                              :: kfl_type  ! Scalar,matrix,vector
     integer(ip)                              :: dim       ! Dimension of properties
     character(5)                             :: name      ! Name of property
     character(5),     pointer                :: wlaws(:)
     integer(ip),      pointer                :: ilaws(:)
     real(rp),         pointer                :: rlaws(:,:)
     type(r1p),        pointer                :: value_ielem(:)
     type(r1p),        pointer                :: value_iboun(:)
     type(r2p),        pointer                :: grval_ielem(:)
     type(r2p),        pointer                :: drval_ielem(:)
     type(r3p),        pointer                :: drval_tur_ielem(:)
     type(r3p),        pointer                :: drval_vel_ielem(:)
     type(r3p),        pointer                :: gdval_ielem(:)
     real(rp),         pointer                :: value_ipoin(:)
     real(rp),         pointer                :: grval_ipoin(:,:)
     real(rp),         pointer                :: value_const(:,:)
     real(rp),         pointer                :: value_default(:)
     type(typ_laws)                           :: llaws(mlaws_ker)
     integer(ip),      pointer                :: update(:,:)
     character (len=5), dimension(:), pointer :: time_function(:)
  end type typ_valpr_ker

  type lis_prope_ker   !Just a list of the properties to improve readability
     type(typ_valpr_ker),pointer :: prop
  end type lis_prope_ker

  !------------------------------------------------------------------------
  !
  ! Physical problem
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_vefun,             &      ! Velocity function
       kfl_tefun,             &      ! Temperature function
       kfl_cofun,             &      ! Concentration function
       kfl_difun,             &      ! Mesh displacement function
       kfl_randseed,          &      ! Enable variable seed for random number generation
       kfl_rough,             &      ! Roughness function: =0: constant,/=0: field
       kfl_canhe,             &      ! Canopy height function: =0: constant,/=0: field
       kfl_heiov,             &      ! Height over terrain function: =0: aprox. by wall distance,/=0: field
       kfl_canla,             &      ! Canopy leaf area density function =0: constant,/=0: field
       kfl_delta,             &      ! Ways to calculate wall distance (=0 prescribed by delta_dom, =1 depending on element size, -1 fix v=0 instead)
       kfl_ustar,             &      ! Law of the wall type
       kfl_walld,             &      ! Wall distance needed
       kfl_walln,             &      ! Wall normal needed
       kfl_walld_field(2),    &      ! Fields assignement to the wallo and wallcoor
       kfl_suppo,             &      ! Support geometry for MM
       kfl_extro,             &      ! Roughness extension
       kfl_prope,             &      ! Properties computed by kermod
       kfl_kxmod_ker,         &      ! k -eps or k-w specific model flag (variable copied to kernel)
       kfl_logva,             &      ! unknowns are logarithmic of usual unknowns (U=log u)
       kfl_noslw_ker,         &      ! Wall law increasing viscosity plus no slip
       kfl_waexl_ker,         &      ! Use of exchange location for wall law
       kfl_waexl_imp_ker,     &      ! Exchange location strategy is performed implicit
       kfl_twola_ker,         &      ! Two layer wall model using auxiliary RANS simulation
       kfl_conta,             &      ! PDN contact
       krestr_2_codno(9),     &      ! Codno to restrict to identify boundaries for two-layer model
       nrestr_2_codno,        &      ! Number of codno ...
       kfl_wlaav_ker,         &      ! Time averaging of ustar for wall law for LES
       kfl_algebra_operations,&      ! Testing of basic algebraic operations
       kfl_vector_size,       &      ! Vector size option
       kfl_temper_vect,       &      ! Activates vectorization in temper
       kfl_aveme_ker                 ! Time averaging use average - average
  !$acc declare create(kfl_ustar)

  integer(ip), pointer ::     &
       kfl_boexc_ker(:)       !(mcodb+1)        ! # Boundary codes at which to calculate exchange point

  real(rp)                 :: &
       denme,                 &      ! Medium density
       visme,                 &      ! Viscosity medium
       gasco,                 &      ! Gas constant R
       conce_relhu,           &      ! Liquid mass fraction at given relative humidity for initializing conce
       u_ref,                 &      ! Reference velocity
       h_ref,                 &      ! Reference height
       k_ref,                 &      ! Reference roughness
       usref,                 &      ! Reference ustar
       windg(3),              &      ! Geostrophic wind (CFDWind2 model)
       delta_dom,             &      ! Wall distance for law of the wall
       delmu_dom,             &      ! Multiplier factor for variable wall distance for law of the wall
       rough_dom,             &      ! Wall distance for law of the wall
       grnor,                 &      ! Gravity norm
       gravi(3),              &      ! Gravity vector
       thicl,                 &      ! Interface thickness
       cmu_st,                &      ! cmu to imposse ABL b.c in nastin
       cpres_ker(2),          &      ! C constant for the windkessel pressure model
       rpres_ker(2),          &      ! R constant for the windkessel pressure model
       dexlo_ker,             &      ! Distance of exchange location method
       tpeav_ker,             &      ! Time period for time-averaging the velocity in wall law
       tlape_ker                     ! Time-averaging period for the two-layer model

  type custo_CFDWind
     integer(ip)  :: kfl_model             ! 1=CFDWind1, 2=CFDWind2
     real(rp)     :: u_ref                 ! Reference velocity
     real(rp)     :: h_ref                 ! Reference height
     real(rp)     :: k_ref                 ! Reference roughness
     real(rp)     :: ustar_ref             ! Reference ustar
     real(rp)     :: l_monin_obukhov       ! Monin-Obukhov mixing length (CFDWind2 model)
     real(rp)     :: u_geostrohpic_wind(3) ! Geostrophic wind (CFDWind2 model)
     real(rp)     :: wind_angle            ! Wind angle
  end type custo_CFDWind

  type typ_subdomain
     integer(ip)  :: kfl_formulation      ! (ALE/Total Lagrangian)
     integer(ip)  :: rotation_number      ! Rotation function number
     character(5) :: rotation_name        ! Rotation function name
     integer(ip)  :: omega_number         ! Omega function number
     character(5) :: omega_name           ! Omega function name
     real(rp)     :: rotation_axis(3)     ! Rotation axis
     real(rp)     :: rotation_center(3)   ! Rotation center
     real(rp)     :: rotation_angle       ! Rotation angle
     real(rp)     :: omega                ! Rotation omega
     integer(ip)  :: translation_number   ! Translation function number
     character(5) :: translation_name     ! Translation name
     real(rp)     :: translation(3)       ! Translation vector
  end type typ_subdomain

  type(typ_subdomain), pointer  :: subdomain(:)

  type(typ_valpr_ker), target   :: densi_ker
  type(typ_valpr_ker), target   :: visco_ker
  type(typ_valpr_ker), target   :: poros_ker
  type(typ_valpr_ker), target   :: condu_ker
  type(typ_valpr_ker), target   :: sphea_ker
  type(typ_valpr_ker), target   :: dummy_ker
  type(typ_valpr_ker), target   :: turmu_ker
  type(typ_valpr_ker), target   :: absor_ker
  type(typ_valpr_ker), target   :: scatt_ker
  type(typ_valpr_ker), target   :: mixin_ker
  type(typ_valpr_ker), target   :: anipo_ker

  !------------------------------------------------------------------------
  !
  ! Numerical treatment
  !
  !------------------------------------------------------------------------

  integer(ip)                   :: &
       kfl_renumbering_npoin,      & ! Renumbering of nodes
       kfl_renumbering_nelem,      & ! Renumbering of elements
       nsfc_renumbering_npoin,     & ! Renumbering of elements
       ndivi,                      & ! Number of mesh divisions
       multiply_with_curvature,    & ! Divide mesh taking into account a curved geometry defined on each element
       kfl_elndi,                  & ! Multiply using element normal method (only for HEX08 and QUA04 elements)
       kfl_mmpar,                  & ! Mesh multiplication parallel strategy, 0=Default, 1=coordinates , 2=New
       kfl_edge_elements,          & ! If edge elements exist
       kfl_rotation_axe,           & ! Axes of rotation of the mesh
       kfl_graph,                  & ! Compute extended graph (needed for RAS precon.)
       kfl_elm_graph,              & ! Compute element graph with halos
       kfl_data_base_array,        & ! Element data base in arrays
       kfl_lface,                  & ! If list of faces is required: LFACG
       kfl_fv_data,                & ! If finite volume data is required
       kfl_grpro,                  & ! Gradient projection (1=closed,0=open)
       kfl_fixsm,                  & ! Default boundary smoothing
       kfl_matrix_grad,            & ! Gradient matrix
       kfl_conbc_ker,              & ! Kernel functions
       kfl_element_bin,            & ! If element bin is needed
       kfl_quali,                  & ! Mesh quality
       kfl_elses,                  & ! Elsest exists
       ielse(nelse),               & ! Elsest int parameters
       kfl_cutel,                  & ! if cut elements exist
       kfl_hangi,                  & ! If hanging nodes exist
       kfl_lapla,                  & ! If hessian matrices are computed
       kfl_defor,                  & ! Mesh deformation
       kfl_coo,                    & ! COO format
       kfl_ell,                    & ! ELL format
       kfl_full_rows,              & ! Full rows graph
       kfl_element_to_csr,         & ! Element to CSR is saved
       kfl_direct_solver,          & ! Direct solver (Alya, Pastix, MUMPS, WSMP,PWSMP #etc.#)
       npoin_mm,                   & ! Support geometry for MM: number of nodes
       nboun_mm,                   & ! Support geometry for MM: number of boundaries
       mnodb_mm,                   & ! Support geometry for MM: max number of nodes per boundary
       number_space_time_function, & ! Number of space/time functions
       number_time_function,       & ! Number of time functions
       number_windk_systems,       & ! Number of Windkessel systems
       number_lookup_tab,          & ! Number of lookup tables
       number_lookup_fw,           & ! Number of lookup frameworks
       deformation_steps,          & ! Number of deformation steps (when KFL_DEFOR /= 0 )
       deformation_strategy,       & ! Deformation strategy
       kfl_duatss,                 & ! Dual time step to precondition transient systems
       fact_duatss,                & ! factor for dual time step
       kfl_conma,                  & ! If consistent mass is needed
       kfl_conma_weighted,         & ! If weighted consistent mass is needed
       kfl_approx_inv_mass,        & ! Approximate inverse mass matrix
       kfl_dmass                     ! Lumped diagonal mass

  real(rp)                 ::      &
       relse(nelse),               & ! Elsest real parameters
       rotation_angle                ! Angle of rotation of the mesh

  integer(ip),    pointer  ::      &
       lnodb_mm(:,:),              & ! Support geometry for MM: connectivity
       ltypb_mm(:),                & ! Support geometry for MM: connectivity
       kfl_type_time_function(:,:)   ! Time function type
  real(rp),       pointer  ::      &
       coord_mm(:,:),              &  ! Support geometry for MM: coordinates
       velav_ker(:,:,:)               ! Time averaged velocity for wall law
  type(bc_nodes), pointer  ::      &
       tncod_ker(:)                  ! Node code type
  type(bc_nodes), pointer  ::      &
       tgcod_ker(:)                  ! Geometrical node code type
  type(bc_bound), pointer  ::      &
       tbcod_ker(:)                  ! Boundary code type
  type(typ_space_time_function), target :: space_time_function(max_space_time_function)  ! Space time functions
  type(typ_time_function),       target :: time_function(max_time_function)              ! Time function
  type(typ_windk_system),        target :: windk_systems(max_windk_systems)              ! Windkessel function
  type(typ_lookup_table),        target :: lookup_tab(max_lookup_tab)                    ! Lookup tables 
  type(typ_tab_coord),           target :: lookup_coords(max_lookup_dim,max_lookup_tab)  ! Coordinates of lookup tables
  type(typ_lookup_framework),    target :: lookup_fw(max_lookup_fw)                      ! Lookup frameworks

  !------------------------------------------------------------------------
  !
  ! Output and postprocess
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       mwitn,                 &      ! Max # witness points
       nwitn,                 &      ! # witness points (local)
       nwitn_all,             &      ! # witness points (global)
       nwitg,                 &      ! # witness geometries (global)
       nwith,                 &      ! # witness meshes (global)
       kfl_posdo,             &      ! Domain postprocess
       kfl_posdi,             &      ! Postprocess division
       kfl_oumes(4),          &      ! Output Mesh needs to be written: geo/set/fixity/sets
       kfl_oustl,             &      ! Output boundary mesh in STL format
       kfl_timeline,          &      ! Output of timeline
       npp_stepo,             &      ! Step over postprocessing frequency
       nfilt,                 &      ! Number of filters
       kfl_abovx,             &      ! Automatic voxel bounding box
       resvx(3),              &      ! Voxel resolution
       kfl_vortx,             &      ! flag of VORTX output
       kfl_vortx_thres,       &      ! flag of threshold VORTX output
       kfl_detection,         &      ! Automatic detection of pattern
       kfl_pixel,             &      ! Output pixel
       plane_pixel,           &      ! Plane along x, y or z
       variable_pixel,        &      ! Variable to postprocess
       number_pixel(2),       &      ! Pixel number
       nsteps_ensemble               ! Number of steps for ensemble (if 0 no ensemble is used)
  real(rp)                 :: &
       bobvx(2,3),            &      ! Voxel Bounding box corners
       travx(3),              &      ! Translation
       thr_veloc,             &      ! threshold VORTX with velocity
       thr_qvort,             &      ! threshold VORTX with qvorticity
       detection_length,      &      ! Characteristic length for detection
       detection_velocity,    &      ! Characteristic velocity for detection
       coord_plane_pixel             ! Plane coordinate for pixel value
  real(rp),    pointer     :: &
       cowit_origi(:,:)              ! Witness points original coordinates
  type(tyfil), target      :: &
       filte(mfilt),          &      ! Filter type
       parfi(10)                     ! Number of parameters
  type(witness_geo), pointer :: &
       gewit(:)                      ! Witness geometry

  type typ_witness_mesh
     type(mesh_type_basic)  :: mesh     ! Witness mesh
     type(witness_geo)      :: geom     ! Witness mesh geometry
     character(200)         :: name     ! Witness mesh name
     type(typ_interp)       :: inte     ! Witness interpolation
  end type typ_witness_mesh
  type(typ_witness_mesh), pointer :: &
       witness_mesh(:)               ! Witness mesh 
  
  !------------------------------------------------------------------------
  !
  ! Optimization and adjoint
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_cost_type,                 &     ! Functional type
       kfl_dvar_type,                 &     ! Design variable type
       kfl_calc_sens,                 &     ! Flag for calculating sensitivity (0 or 1)
       kfl_adj_prob,                 &      ! Adjoint sansitivities (1 = ON)
       kfl_cos_opt,                 &       ! Calculate functional  (1 = ON)
       kfl_nwall,                     &     ! Number of wall nodes
       kfl_ndvars_opt                       ! Number of design variables
  
  real(rp)               :: &
       costf                              ! cost function value (1 shot)

  real(rp), pointer      :: &
       sens_mesh(:,:),                   &      ! mesh sensitivity values respecto to the mesh node coordinates
!        untur_fix(:),                   &      !
       force_norm(:),                   &        ! drag, lift,... forces (normalized)
       sens(:)                              ! sensitivity values

  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       lastm_ker,             &      ! Last solved module
       nvort,                 &      ! number of vortex
       nvort_total,           &      ! number total of vortex
       npoin_total_filt,      &      ! number total of filtered points
       npoin_filt,            &      ! number of filtered points
       nelem_total_filt,      &      ! number total of filtered elements
       nelem_filt,            &      ! number of filtered points
       number_event                  ! Event number
  real(rp)                 :: &
       avwei_ker                     ! weight for time averaging

  real(rp),   pointer      :: &
       cowit(:,:),            &      ! Witness points coordinates
       uwall_ker(:),          &      ! Save wall distance solution
       uwal2_ker(:),          &      ! Save wall distance solution
       shwit(:,:),            &      ! Witness point shape functions
       dewit(:,:,:),          &      ! Witness point shape functions derivatives
       displ_ker(:,:)                ! Displacement support boundary for MM
  integer(ip), pointer     :: &
       lewit(:)                      ! List of witness point element

  !------------------------------------------------------------------------
  !
  ! Matrices
  !
  !------------------------------------------------------------------------

  real(rp),   pointer      ::    &
       matrix_grad(:,:)              ! Gradient matrix

  !------------------------------------------------------------------------
  !
  ! Boundary conditions of wall distance and wall normal problem
  !
  !------------------------------------------------------------------------

  integer(ip), parameter   ::    &
       mfunc_walld_ker=10            ! Maximum number of functions
  integer(ip), pointer     ::    &
       kfl_funno_walld_ker(:),   &   ! Nodal function
       kfl_funbo_walld_ker(:),   &   ! Boundary function
       kfl_fixno_walld_ker(:,:), &   ! Nodal fixity
       kfl_fixbo_walld_ker(:),   &   ! Boundary fixity
       kfl_funty_walld_ker(:,:), &   ! Function type and number of parameters
       kfl_fixno_walln_ker(:,:), &   ! Nodal fixity
       kfl_fixbo_walln_ker(:),   &   ! Boundary fixity
       lexlo_ker(:,:)                ! List for the exchange location for wall law
  type(r1p),   pointer     ::    &
       funpa_walld_ker(:)            ! Function parameters
  real(rp),    pointer     ::    &
       bvess_walld_ker(:,:),     &   ! Nodal value
       bvnat_walld_ker(:),       &   ! Boundary value
       bvess_walln_ker(:,:),     &   ! Nodal value
       velel_ker(:,:),           &   ! Velocity at the exchange location for wall law
       temel_ker(:,:),           &   ! Temperature at the exchange location for wall law
       shape_waexl(:,:)              ! shape functions associated to excange location (for implicit)

  type(typ_color_coupling) :: wallcoupling_waexl
  type(typ_color_coupling) :: wallcoupling_extr_boun ! used to extrapolate shear stress from boundaries for two layer model
  integer(ip), pointer     ::    &
       ptb_to_use(:),            &   ! boundary point to use
       is_interior(:)                !
  !------------------------------------------------------------------------
  !
  ! Boundary conditions of deformation problem
  !
  !------------------------------------------------------------------------

  integer(ip), pointer     ::    &
       kfl_funno_defor_ker(:),   &   ! Nodal function
       kfl_fixno_defor_ker(:,:), &   ! Nodal fixity
       kfl_funty_defor_ker(:,:)      ! Function type and number of parameters
  real(rp),    pointer     ::    &
       bvess_defor_ker(:,:)          ! Nodal value

  !------------------------------------------------------------------------
  !
  ! Boundary conditions of roughness extension problem
  !
  !------------------------------------------------------------------------

  integer(ip), parameter   ::    &
       mfunc_rough_ker=10            ! Maximum number of functions
  integer(ip), pointer     ::    &
       kfl_funno_rough_ker(:),   &   ! Nodal function
       kfl_funbo_rough_ker(:),   &   ! Boundary function
       kfl_fixno_rough_ker(:,:), &   ! Nodal fixity
       kfl_fixbo_rough_ker(:),   &   ! Boundary fixity
       kfl_funty_rough_ker(:,:)      ! Function type and number of parameters
  type(r1p),   pointer     ::    &
       funpa_rough_ker(:)            ! Function parameters
  real(rp),    pointer     ::    &
       bvess_rough_ker(:,:),     &   ! Nodal value
       bvnat_rough_ker(:)            ! Boundary value

  !------------------------------------------------------------------------
  !
  ! Boundary conditions of support geometry
  !
  !------------------------------------------------------------------------

  integer(ip), pointer     ::    &
       kfl_fixno_suppo_ker(:,:)      ! Nodal fixity
  real(rp),    pointer     ::    &
       bvess_suppo_ker(:,:)          ! Nodal value

  !------------------------------------------------------------------------
  !
  ! Aditional scalars & arrays for no slip wall law
  !
  !------------------------------------------------------------------------

  integer(ip)              ::    &
       kount_nsw_ele_ker,        &       ! # elements that intervine in no slip wall law
       kount_nsw_poi_ker                 ! # nodes    that intervine in no slip wall law

  integer(ip), pointer     ::    &
       kfl_nswel_ker(:),         &       ! list of elements that intervine in no slip wall law
       kfl_nswpo_ker(:)                  ! list of node     that intervine in no slip wall law

  real(rp),    pointer     ::    &
       normal_nsw_ker(:,:),      &       ! element normal
       avta1_nsw_ker(:,:),       &       ! element normal
       fact_nsw_ker(:)                   ! correction factor for the difference between the variational traction and the one obtained from vel grad

  integer(ip), pointer     ::    &
       kfl_fixbo_nsw_ker(:)              ! Boundary fixity

  real(rp),    pointer     :: &
       avupo_ker(:,:)                    ! Post-processing time averaged velocity wall law  - used not only for no slip wall
  !
  ! All this part might deserve better names and actually it is may not be necesaryy to use a derived data type
  !
  type exch_loc_elemv
     integer(ip)           :: nbogp
     real(rp)              :: fact
     real(rp)              :: vel_aux(3)
     real(rp)              :: velav(3)
  end type exch_loc_elemv
  type(exch_loc_elemv),   pointer  :: &
       lnsw_exch(:)                  ! For obtaining average(over all boundary gp) exchange location average (over time) velocity related to gauss points corresponding to an element

  !------------------------------------------------------------------------
  !
  ! Events
  !
  !------------------------------------------------------------------------

  character(50) :: events_directory

  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------

  real(rp),    pointer     ::    &
       el_nsw_visc(:)        ! elemental aditional viscosity for no slip wall law

end module def_kermod
