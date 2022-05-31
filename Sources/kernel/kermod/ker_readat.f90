!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_readat.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Read kermod data.
!> @details Read kermod data.
!> @}
!-----------------------------------------------------------------------

subroutine ker_readat()
  use def_kintyp
  use def_elmtyp
  use def_master
  use def_kermod
  use def_inpout
  use def_domain
  use mod_witness
  use def_solver,                  only : SOL_DIRECT_SOLVER_MUMPS
  use def_solver,                  only : SOL_DIRECT_SOLVER_ALYA
  use def_solver,                  only : SOL_DIRECT_SOLVER_SKYLINE_ALYA
  use def_solver,                  only : SOL_DIRECT_SOLVER_PASTIX
  use def_solver,                  only : SOL_DIRECT_SOLVER_GPUQR
  use def_solver,                  only : SOL_DIRECT_SOLVER_WSMP
  use def_solver,                  only : SOL_DIRECT_SOLVER_PWSMP
  use mod_ker_space_time_function, only : space_time_function_number
  use mod_ker_polynomial,          only : ker_polynomial_readat, ker_polynomial_name
  use mod_opebcs,                  only : opebcs_initialization_structure
  use mod_opebcs,                  only : opebcs_initialization_variable
  use mod_ecoute,                  only : ecoute
  use mod_ker_regularization,      only : kfl_regularization, set_regularization_type, reg_exp
  use mod_ker_regularization,      only : reg_exp_linear, reg_exp_quadratic, reg_exp_garantzha, reg_exp_traslation
  use mod_ker_regularization,      only : reg_exp_compression, reg_exp_quadratic_comp, reg_exp_log, reg_exp_log
  use mod_ker_regularization,      only : reg_exp_cub, kfl_second
  use mod_ker_tendencies,          only : kfl_tendencies_ker, read_tendencies
  use mod_ker_subdomain,           only : ker_subdomain_read_data
  use mod_ker_subdomain,           only : ker_subdomain_function_name_to_number
  use mod_messages,                only : livinf
  use mod_messages,                only : messages_live 
  use def_master,                  only : intost
  use def_master,                  only : retost
  use mod_memory,                  only : memory_alloca
  use mod_read_domain_arrays,      only : read_domain_arrays_types
  use mod_physics,                 only : physics_humidity
  use mod_ker_proper,              only : MATRIX_PROPERTY
  use mod_domain,                  only : domain_memory_allocate
  use mod_witness,                 only : witness_geometry_initialization
  use mod_output_postprocess,      only : output_postprocess_read_geometry
  use mod_ker_discrete_function,   only : ker_discrete_function_read
  use mod_AMR,                     only : AMR_readat
  use mod_interp_tab,              only : tab_load_file,fw_allocate
  implicit none

  integer(ip)  :: imodu,jmodu,idime,imate,ipara,ipoin,ifunc,kdime,iauxi,ii,jj
  integer(ip)  :: iboun,dummi,ktype,lexis_mm(nelty),ifunp,nfunp,iwitg
  integer(ip)  :: nposi,icosi,order_loc
  integer(ip)  :: i_lookup_fw, i_lookup_tab 
  real(rp)     :: dummr,scawi(3),tem_relhu,relhu
  character(5) :: icha1='     ', icha2='     '
  character(5) :: wfname
  character(5) :: wfnamt
  character(5) :: wfnamc
  character(5) :: wfnamd
  character(5) :: wcsys
  type(typ_lookup_table),     pointer :: ptr_lookup_tab                    
  type(typ_tab_coord),        pointer :: ptr_lookup_coords(:)  
  type(typ_lookup_framework), pointer :: ptr_lookup_fw

  call memory_alloca(mem_modul(1:2,modul), 'KFL_BOEXC_KER', 'ker_memory', kfl_boexc_ker, mcodb+1)

  !-----------------------------------------------------------------------
  !
  ! Boundary conditions: allocate and initialize
  !
  ! 1. Mesh deformation/smoothing
  ! 2. Wall distance
  ! 3. Deformation
  ! 4. Mesh support surface
  ! 5. Wall normal
  !
  !-----------------------------------------------------------------------
  !
  ! Node
  !
  call opebcs_initialization_structure(7_ip,tncod_ker)
  call opebcs_initialization_variable (1_ip,tncod_ker)
  !
  ! Boundary
  !
  call opebcs_initialization_structure(7_ip,tbcod_ker)
  call opebcs_initialization_variable (1_ip,tbcod_ker)
  !
  ! Geometrical
  !
  call opebcs_initialization_structure(2_ip,tgcod_ker)
  call opebcs_initialization_variable (1_ip,tgcod_ker)

  if( INOTSLAVE ) then
     !
     ! Physical problem
     !
     kfl_vefun     =  0                                     ! Velocity function
     kfl_tefun     =  0                                     ! Temperature function
     kfl_cofun     =  0                                     ! Concentration function
     kfl_difun     =  0                                     ! Mesh displacement function
     kfl_randseed  =  1                                     ! Variable seed for random numbers: enabled by default
     kfl_rough     = -1                                     ! No roughness
     kfl_canhe     = -1                                     ! No canopy
     kfl_heiov     = -1                                     ! No height over terrain
     kfl_canla     =  0                                     ! Constant leaf area density (canopy)
     kfl_delta     =  0                                     ! Ways to calculate wall distance
     kfl_ustar     =  0                                     ! Ways to calculate U* (=0: Reichart, 1=ABL, 2=ABL2)
     kfl_walld     =  0                                     ! Wall distance not needed
     kfl_walln     =  0                                     ! Wall normal not needed
     kfl_suppo     =  0                                     ! Support geometry for MM
     kfl_extro     =  0                                     ! Roughness extension
     kfl_prope     =  0                                     ! Kermod in charge of properties
     kfl_kxmod_ker =  0                                     ! kermod k eps model modification
     kfl_noslw_ker =  0                                     ! Do not use wall law increasing viscosity plus no slip
     kfl_waexl_ker =  0                                     ! Do not use exchange location for wall law
     kfl_waexl_imp_ker   = 0                                ! Exchange location treated explicit by default
     kfl_twola_ker       = 0                                ! No auxiliary RANS simulation used
     krestr_2_codno(1:9) = 0                                ! codno to restrict to identify boundaries for two-layer model
     nrestr_2_codno      = 0                                ! number of codno ...
     kfl_wlaav_ker       = 0                                ! Do not use time-averaging for wall law
     kfl_aveme_ker       = 1                                !
     kfl_conta     = 0                                      ! PDN contact
     
     denme         = -1.0_rp                                ! Medium density
     visme         = -1.0_rp                                ! Viscosity medium
     gasco         =  8.3144621_rp                          ! Gas constant
     u_ref         = -1.0_rp                                ! Reference velocity
     h_ref         = -1.0_rp                                ! Reference height
     k_ref         = -1.0_rp                                ! Reference roughness
     usref         = -1.0_rp                                ! Reference ustar
     windg         =  0.0_rp                                ! Geostrophic wind (CFDWind2 model)
     delta_dom     =  0.0_rp                                ! Distance to the wall
     delmu_dom     =  1.0_rp                                ! Multiplier factor for variable wall distance for law of the wall
     rough_dom     =  0.0_rp                                ! Wall roughness
     grnor         =  0.0_rp                                ! Gravity norm
     gravi         =  0.0_rp                                ! Gravity vector
     thicl         = -1.0_rp                                ! Thickness of free surface
     cmu_st        =  0.0_rp                                ! Cmu for turbul k-eps, needed in Nastin b.c
     dexlo_ker     =  0.0_rp                                ! Exchange location distance
     tpeav_ker     =  0.0_rp                                ! Period for time-averaging the velocity in wall law
     tlape_ker     =  0.0_rp                                ! Time-averaging period for two-layer model
     kfl_algebra_operations = 2                             ! Testing of basic algebraic operations
     cpres_ker     =  0.0_rp                                ! C constant for the windkessel pressure model
     rpres_ker     =  0.0_rp                                ! R constant for the windkessel pressure model
     num_lobas     =  0                                     ! Nu,ber of coordinate systems

     kfl_boexc_ker = 1                                      ! With this it calculates exchange points for all boundaries ( Matias was forced to do this when he moved to kernel)
     ! Now we will be able to tell Alya to do it only for certain boundaries
     !
     ! Numerical treatment
     !
     if( ISEQUEN ) then
        kfl_renumbering_npoin = 0                           ! Renumbering strategy
        kfl_renumbering_nelem = 0                           ! Renumbering strategy
     else
        kfl_renumbering_npoin = 1                           ! Renumbering strategy
        kfl_renumbering_nelem = 1                           ! Renumbering strategy
     end if
     nsfc_renumbering_npoin = 256                           ! Renumbering strategy: default SFC size
     ndivi         = 0                                      ! Number of divisions
     multiply_with_curvature = 0                            ! Multiply with curvature or not
     kfl_elndi     = 0                                      ! Multiply using element normal method (only for HEX08 and QUA04)
     kfl_mmpar     = 0                                      ! Mehs mutliplication strategy
     kfl_edge_elements = 0                                  ! No edge elements
     kfl_rotation_axe = 0                                   ! Axes of rotation of the mesh
     kfl_graph     = 0                                      ! Extended graph
     kfl_elm_graph = 0                                      ! Element graph with halos
     kfl_savda     = 0                                      ! Do not save element data base
     kfl_data_base_array = 0                                ! Save elemental element functions in an array (used for GPU)
#ifdef VECTOR_SIZE_VARIABLE
     VECTOR_SIZE   = 1
#endif
     kfl_vector_size = 0                                    ! Vector size option
     kfl_temper_vect = 0_ip                                    ! Vector size option
     kfl_lface     = 0                                      ! If list of faces is required: LFACG
     kfl_fv_data   = 0                                      ! If finite volume data is required
     kfl_grpro     = 1                                      ! Gradient projection (1=closed,0=open)
     kfl_matrix_grad = 0                                    ! Gradient matrix
     kfl_conbc_ker = 1                                      ! Constant kernel b.c.
     kfl_element_bin = 0                                    ! No element bin
     kfl_elses     = 0                                      ! Elsest
     ielse(1)      = 100                                    ! nx
     ielse(2)      = 100                                    ! ny
     ielse(3)      = 100                                    ! nz
     ielse(4)      = 1                                      ! data format (0=type,1=linked list)
     ielse(5)      = 10                                     ! Maximum number of possible meshes
     ielse(6)      = 1                                      ! Not used
     ielse(7)      = 0                                      ! Output unit
     ielse(8)      = 1                                      ! Search strategy (0=bin,1=Quad)
     ielse(9)      = 100                                    ! Points per node for Quad/Octtree
     ielse(10)     = 0                                      ! Result output frequency
     ielse(11)     = 1                                      ! Neighboring-boxes-search radius, 0: search until all boxes finished
     ielse(12)     = 0                                      ! Postprocess mesh
     ielse(13)     = 0                                      ! Postprocess results
     ielse(14)     = 0                                      ! Dont check
     ielse(15)     = 0                                      ! Dont force
     ielse(16)     = 0                                      ! Save element bounding box

     relse(1)      = 1.0e-6_rp                              ! Tolerance for iteration
     rotation_angle = 0.0_rp                                ! Angle of rotation of the mesh (in degrees)

     kfl_cutel     = 0                                      ! If cut elements exist
     kfl_hangi     = 0                                      ! If hanging node exist
     kfl_lapla     = 1                                      ! Laplacians are computed
     kfl_defor     = 0                                      ! Do not deform mesh
     kfl_coo       = 0                                      ! COO format is off
     kfl_full_rows = 0                                      ! Full rows graph is off
     kfl_element_to_csr = 0                                 ! Do not save element to CSR array
     kfl_direct_solver  = SOL_DIRECT_SOLVER_ALYA            ! Direct solver option
     element_bin_boxes(1:3) = 4                             ! Element bin
     npoin_mm      = 0                                      ! Support geometry for MM: number of nodes
     nboun_mm      = 0                                      ! Support geometry for MM: number of boundaries
     mnodb_mm      = ndime                                  ! Support geometry for MM: max number of nodes per boundary
     number_space_time_function = 0                         ! # of space/time functions
     number_time_function       = 0                         ! # of time functions
     number_windk_systems       = 0                         ! # of Windkesssel functions
     number_lookup_tab          = 0                         ! # of lookup tables 
     number_lookup_fw           = 0                         ! # of lookup frameworks 

     do ifunc = 1,max_space_time_function
        space_time_function(ifunc) % ndime = 1              ! Space/time function dimension
        space_time_function(ifunc) % nexpr = 0              ! Space/time function size of the expression
        space_time_function(ifunc) % name  = ' '            ! Space/time function name
     end do
     do ifunc = 1,max_time_function
        time_function(ifunc) % npara = 1                    ! Time function dimension
        time_function(ifunc) % kfl_type = 0                 ! Time function size of the expression
        time_function(ifunc) % name  = ' '                  ! Time function name
     end do
     do ifunc = 1,max_windk_systems
        windk_systems(ifunc) % sysid            = 0         ! ID of the system
        windk_systems(ifunc) % name             = ' '       ! Time function name
        windk_systems(ifunc) % wdks_model       = 0         ! Time function size of the expression
        windk_systems(ifunc) % nparam           = 1         ! Time function dimension
        windk_systems(ifunc) % x_in             = 0.0_rp    ! read value
        windk_systems(ifunc) % y_out            = 0.0_rp    ! computed value
        windk_systems(ifunc) % stored_time_step = 0.0_rp    ! computed value
        windk_systems(ifunc) % iflow_nsi        = 0         ! nastin internal parameter
     end do
     deformation_steps    = 1                               ! Number of deformation steps (when KFL_DEFOR /= 0 )
     deformation_strategy = 6                               ! Deformation strategy
     kfl_duatss = 0                                         ! No dual time stepping as preconditioner
     fact_duatss = 10                                       ! Make time step 10 times smaller
     kfl_conma            = 0                               ! Consistent mass matrix
     kfl_conma_weighted   = 0                               ! Weighted Consistent mass matrix
     kfl_approx_inv_mass  = 0                               ! Approx invers
     kfl_dmass            = 0                               ! Diagonaly mass matrix
     kfl_logva = 0                                          ! logarithmic variable
     kfl_regularization = .false.
     kfl_second = .false.
     kfl_tendencies_ker =.false.                            ! Meso-to-micro coupling

     !
     ! Optimization
     !
     kfl_cost_type  = 0                                     ! Functional type
     kfl_dvar_type  = 0                                     ! Design variable type
     kfl_adj_prob   = 0                                     ! Adjoint sansitivities (1 = ON)
     kfl_ndvars_opt = 0                                     ! Number of design variables
     kfl_cos_opt    = 0
     !
     ! Output and postprocess
     !
     nwitn_all          =  0                                ! # witness points (read by master)
     nwitn              =  0                                ! # witness points
     mwitn              =  50                               ! Max # witness points
     nwitg              =  0                                ! Number of witness geometries
     nwith              =  0                                ! Number of witness meshes
     kfl_posdo          =  1                                ! Domain postprocess (default)
     kfl_posdi          =  0                                ! Postprocess on original level
     kfl_oumes          =  0                                ! Do not output mesh
     kfl_oumes(1)       =  1                                ! Output mesh: only geometry
     kfl_oustl          =  0                                ! Output boundary mesh STL format
     kfl_quali          =  0                                ! Mesh quality
     npp_stepo          = -1                                ! Do not step over the defined postprocesing steps
     nfilt              =  0                                ! # filters
     kfl_abovx          =  0                                ! Automatic voxel bounding box
     resvx              =  32                               ! Voxel resolution: 32 x 32 x 32
     kfl_vortx          =  0                                ! Vortex extraction option off
     kfl_vortx_thres    =  0                                ! Vortex threshold option off
     kfl_detection      =  0                                ! Automatic detection
     kfl_livee          =  1                                ! Live info each [kfl_livee] steps
     kfl_pixel          =  0                                ! Output ppm format
     plane_pixel        =  1                                ! Plane along x,y,z
     variable_pixel     =  0                                ! Variable to postprocess
     number_pixel       =  50                               ! Numebr of pixels
     nsteps_ensemble    =  0                                ! No ensemble is used
     thr_veloc          =  0.1_rp                           ! Vortex threshold velocity
     thr_qvort          =  10.0_rp                          ! Vortex threshold qvort
     travx              =  0.0_rp                           ! Translation
     bobvx              =  0.0_rp                           ! Voxel Bounding box corners
     scawi              =  1.0_rp                           ! Scale for witness points
     detection_length   =  1.0_rp                           ! Characteristic length for detection
     detection_velocity =  1.0_rp                           ! Characteristic velocity for detection
     coord_plane_pixel  =  0.0_rp                           ! Coordinate of the plane

     !-------------------------------------------------------------------
     !
     ! Read/write unit
     !
     !-------------------------------------------------------------------

     if( momod(modul) % lun_pdata <= 0 ) then
        call livinf(0_ip,'KERMOD FILE DOES NOT EXITS: USE DEFAULT OPTIONS',0_ip)
        return
     end if

     lispa = 0
     lisda = momod(modul) % lun_pdata ! Reading file
     lisre = momod(modul) % lun_outpu ! Writing file
     !--><group>
     !-->       <groupName>Physical_problem</groupName>
     !-->       <subGroup>
     !-->           <inputLine>
     !-->               <inputLineName>PRE_PENALIZATION</inputLineName>
     !-->               <inputLineHelp>Penalize pressure in time. To be used when converging to a steady state when the pressure presents an odd-even decoupling in time. This option helps converging by adding a damping effect.</inputLineHelp>
     !-->               <inputElement>
     !-->                   <inputElementType>edit</inputElementType>
     !-->                   <inputLineEditName>VALUE</inputLineEditName>
     !-->                   <inputLineEditValue>5</inputLineEditValue>
     !-->               </inputElement>
     !-->           </inputLine>
     !-->       </subGroup>
     !-->   </group>
     !
     !.md<module>kernel
     !.md<input>case.ker.dat
     !.md<pos>0
     !.md<sec>
     !.md# Physical problem
     !.mdDefine the physical problem of kermod module. This field contains some variables shared by
     !.mdthe modules, like physical properties, distance to the wall, wall laws, etc.
     !.md<code>
     !.md<0>PHYSICAL_PROBLEM
     !.md<com>
     !
     call ecoute('ker_readat')
     do while( words(1) /= 'PHYSI' )
        call ecoute('ker_readat')
     end do
     call ecoute('ker_readat')

     do while(words(1) /= 'ENDPH' )

        if(      words(1) == 'SUBDO' ) then
           !
           !.md<1>SUBDOMAIN:
           !.md<field>SUBDOMAIN
           !.md<com>_..._
           !
           call ker_subdomain_read_data()

        else if( words(1) == 'VELOC' ) then
           !
           !.md<1>VELOCITY: OFF/HELIC/FUNCTION/SPACE_TIME_FUNCTION/FIELD= int
           !.md<field>VELOCITY
           !.md<com>In case a velocity field is needed, for example to transport particles
           !.md<com>with PARTIS module, and no module solves for the velocity, one can define a velocity
           !.md<com>field.
           !
           if(      words(2) == 'OFF  ' ) then
              kfl_vefun = 0
           else if( words(2) == 'NULL ' .or. words(2) == 'ZERO ' ) then
              kfl_vefun = 3
           else if( words(2) == 'HELIC' ) then
              kfl_vefun = 7
           else if( words(2) == 'FIELD' ) then
              kfl_vefun = -getint('FIELD',-1_ip,'#VELOCITY FROM A FIELD')
           else if( exists('SPACE') ) then                                           ! >1000: Space time function
              wfname = getcha('SPACE','NULL ','#Space/time Function name')
              kfl_vefun = 1001
           else
              if( exists('FUNCT') ) &
                   kfl_vefun = getint('FUNCT',1_ip,'#VELOCITY FUNCTION FOR LAGRANGIAN PARTICLE')
           end if

        else if( words(1) == 'DISPL' ) then
           !
           !.md<1>DISPLACEMENT: OFF/SPACE_TIME_FUNCTION/FIELD= int
           !.md<field>DISPLACEMENT
           !.md<com>In case a displcaement field is needed, for example to transport particles
           !.md<com>with NASTIN module, and no module solves for the mesh displacement, one can define a displacement
           !.md<com>field.
           !
           if(      words(2) == 'OFF  ' ) then
              kfl_difun = 0
           else if( words(2) == 'NULL ' .or. words(2) == 'ZERO ' ) then
              kfl_difun = 3
           else if( words(2) == 'FIELD' ) then
              kfl_difun = -getint('FIELD',-1_ip,'#VELOCITY FROM A FIELD')
           else if( exists('SPACE') ) then                                           ! >1000: Space time function
              wfnamd = getcha('SPACE','NULL ','#Space/time Function name')
              kfl_difun = 1001
           else
              if( exists('FUNCT') ) &
                   kfl_difun = getint('FUNCT',1_ip,'#VELOCITY FUNCTION FOR LAGRANGIAN PARTICLE')
           end if

        else if( words(1) == 'TEMPE' ) then
           !
           ! ADOC[1]> TEMPER: FIELD=
           !
           if( words(2) == 'FIELD' ) then
              kfl_tefun = -getint('FIELD',-1_ip,'#TEMPERATURE FROM A FIELD')
           else if( words(2) == 'NULL ' .or. words(2) == 'ZERO ' ) then
              kfl_tefun = 3
           else if( exists('SPACE') ) then                                           ! >1000: Space time function
              wfnamt = getcha('SPACE','NULL ','#Space/time Function name')
              kfl_tefun = 1001
           else
              if( exists('FUNCT') ) &
                   kfl_tefun = getint('FUNCT',1_ip,'#TEMPERATURE FUNCTION FOR LAGRANGIAN PARTICLE')
           end if

        else if( words(1) == 'CONCE' ) then
           !
           ! ADOC[1]> CONCE: FIELD=
           !
           if( words(2) == 'FIELD' ) then
              kfl_cofun = -getint('FIELD',-1_ip,'#CONCENTRATION FROM A FIELD')
           else if( words(2) == 'NULL ' .or. words(2) == 'ZERO ' ) then
              kfl_cofun = 3
           else if( exists('SPACE') ) then                                           ! >1000: Space time function
              wfnamc = getcha('SPACE','NULL ','#Space/time Function name')
              kfl_cofun = 1001
           else
              if( exists('FUNCT') ) &
                   kfl_cofun = getint('FUNCT',1_ip,'#CONCENTRATION FUNCTION FOR LAGRANGIAN PARTICLE')

              if( exists('RELHU') ) then
                  !
                  ! Water concentration based on relative humidity
                  !
                  kfl_cofun     = 666_ip
                  tem_relhu     = getrea('TEMPE',298.15_rp,'#DEFAULT TEMPERATURE TO EVALUATE RELATIVE HUMIDITY')
                  relhu         = getrea('RELHU',0.0_rp,   '#DEFAULT RELATIVE HUMIDITY')
                  conce_relhu   = physics_humidity(tem_relhu,relhu)

              endif
           end if

       else if( words(1) == 'RANDO' ) then
           !
           ! Random seed treatment
           !
           if( (words(2) == 'OFF  ') .or. (words(2) == 'NONE ') .or. (words(2) == 'DISAB') ) then
               kfl_randseed = 0
           end if

        else if( words(1) == 'PROPE' ) then
           !
           !.md<1>PROPERTIES
           !.md<2>MATERIAL= int
           !.md<3>DENSITY:       CONSTANT/BIFLUID/LOW_MACH/K_LOW_MACH,MIXTURE  PARAMETERS= real1, real2..., DEFAULT_VALUE= real3
           !.md<3>VISCOSITY:     CONSTANT/BIFLUID/SUTHERLAND/MU_MIXTURE,       PARAMETERS= real1, real2..., DEFAULT_VALUE= real3
           !.md<3>POROSITY:      CONSTANT/CANOPY/DAFOR/TIME_FUNCTION,          PARAMETERS= real1, real2..., DEFAULT_VALUE= real3, TIME_FUNCTION = char1
           !.md<3>CONDUCTIVITY:  CONSTANT/SUTHERLAND/K_MIXTURE,                PARAMETERS= real1, real2..., DEFAULT_VALUE= real3
           !.md<3>SPECIFIC_HEAT: CONSTANT/CP_MIXTURE,                          PARAMETERS= real1, real2..., DEFAULT_VALUE= real3
           !.md<3>TURBUL       : WALE/STDKE/STTKO/KOMEG,    UPDP1 = ENDITE/ENDINN/ENDSTE/DOITER/ENDBLCK/NEVER  UPDP2 = ENDITE/ENDINN/ENDSTE/DOITER/ENDBLCK/NEVER
           !.md<2>END_MATERIAL
           !.md<1>END_PROPERTIES
           !.md<field>PROPERTIES
           !.md<com>For each material (MATERIAL= int), the properties laws and values are defined.
           !.md<com>If the property is constant, then only one parameter is required.
           !.md<com>For example: DENSITY: CONSTANT, PARAMETERS=1.2047
           !.md<com>Default values may be required because properties can be needed before
           !.md<com>for computing the initial values of the unknowns
           !
           ! Properties
           !
           call ecoute('ker_readat')
           imate = 1
           do while( words(1) /= 'ENDPR' )
              if( words(1) == 'MATER' ) then
                 imate = getint('MATER',1_ip,'#MATERIAL NUMBER')
                 if( imate < 1 .or. imate > nmate ) call runend('KER_READAT: WRONG MATERIAL')
              end if

              if( words(1) == 'DENSI' ) then
                 !
                 ! Density
                 !
                 densi_ker % kfl_exist                = 1
                 if (exists('NASTA')) then
                    densi_ker % wlaws(imate)          = 'DNAST'
                 else if (exists('LOWMA') .and. exists('MIXTU')) then
                    densi_ker % wlaws(imate)          = 'KLOWM'
                 else if (exists('LOWMA') .and. exists('GAUSS')) then
                    densi_ker % wlaws(imate)          = 'LOWMG'
                 else if (exists('LOWMA') .and. exists('GMIXT')) then
                    densi_ker % wlaws(imate)          = 'GKLOW'
                 else if (exists('LOWMA') .and. exists('GPLOO')) then
                    densi_ker % wlaws(imate)          = 'LMGPL'
                    !
                 else if (exists('SPRAY')) then
                    densi_ker % wlaws(imate)          = 'SPRAY'
                 else if (exists('CMCGP')) then
                    densi_ker % wlaws(imate)          = 'CMCGP'
                 else if (exists(ker_polynomial_name)) then
                    call ker_polynomial_readat( densi_ker, imate )
                    !
                 else
                    densi_ker % wlaws(imate)          = words(2)
                 endif
                 densi_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 densi_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT DENSITY')

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, densi_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, densi_ker%update(2,imate))
                 end if

              else if( words(1) == 'VISCO' ) then
                 !
                 ! Viscosity
                 !
                 visco_ker % kfl_exist                = 1
                 if (exists('MIXTU')) then
                    visco_ker % wlaws(imate)          = 'MUMIX'
                 else if (exists('MUGPT')) then
                    visco_ker % wlaws(imate)          = 'MUGPT'
                 else if (exists('TABLE')) then
                    visco_ker % wlaws(imate)          = 'MUTAB'
                 else if (exists('SUTHE') .and. exists('GAUSS')) then  ! Modifier GAUSS_POINT goes at end of line
                    visco_ker % wlaws(imate)          = 'GPSUT'
                 else if (exists('ABL  ') ) then                       ! Atmospheric boundary layer
                    visco_ker % wlaws(imate)          = 'ABL  ' 
                 else if (exists('SPRAY') ) then                       ! Spray model: liquid/gas multiphase flow
                    visco_ker % wlaws(imate)          = 'SPRAY'
                 else if (exists('GAUSS')) then                        ! Gaussian function to generate fluctuations by digital filters
                    visco_ker % wlaws(imate)          = 'GAUSS'        ! with diffusion
                    !
                 else if (exists('NAHME')) then                        ! Gaussian function to generate fluctuations by digital filters
                    visco_ker % wlaws(imate)          = 'NAHME'        ! Nahme's approximation
                 else if (exists(ker_polynomial_name)) then
                    call ker_polynomial_readat( visco_ker, imate )
                    !
                 else
                    visco_ker % wlaws(imate)          = words(2)
                 end if
                 visco_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 visco_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT VISCOSITY')

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, visco_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, visco_ker%update(2,imate))
                 end if

              else if( words(1) == 'POROS' ) then
                 !
                 ! Porosity
                 !
                 poros_ker % kfl_exist                = 1
                 poros_ker % wlaws(imate)             = words(2) ! CONST | CANOPY | FUNCTION | VALVE
                 poros_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 poros_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT POROSITY')
                 poros_ker % time_function(imate)     = getcha('TIMEF','NULL ','#DEFAULT TIME FUNCTION')

                 !check parameters for canopy (Height, CD, LAD, Zmax/h )

                 if ( poros_ker % wlaws(imate)=='CANOP'.and. &
                      (poros_ker % rlaws(4,imate) < -epsilon(1.0_rp).or.&
                      poros_ker%rlaws(4,imate)    > 1.0_rp+epsilon(1.0_rp)) ) &
                      call runend('KER_READAT:WRONG CANOPY ZMAX PARAM, SHOULD BE BETWEEN 0.0 AND 1.0' )

                 ! checks when to update
                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, poros_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, poros_ker%update(2,imate))
                 end if

              else if( words(1) == 'CONDU' ) then
                 !
                 ! Conductivity
                 !
                 condu_ker % kfl_exist                = 1
                 if (exists('MIXTU')) then
                    condu_ker % wlaws(imate)          = 'KMIXT'
                 else if (exists('KGPTA')) then
                    condu_ker % wlaws(imate)          = 'KGPTA'
                 else if (exists('TABLE')) then
                    condu_ker % wlaws(imate)          = 'KTABL'
                 else if (exists('SUTHE') .and. exists('GAUSS')) then  ! Modifier GAUSS_POINT goes at end of line
                    condu_ker % wlaws(imate)          = 'GPSUT'
                 else if (exists('CMCNP')) then
                    condu_ker % wlaws(imate)          = 'CMCNP'
                    !
                 else if (exists(ker_polynomial_name)) then
                    call ker_polynomial_readat( condu_ker, imate )
                    !
                 else
                    condu_ker % wlaws(imate)          = words(2)
                 end if
                 condu_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 condu_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT CONDUCTIVITY')

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, condu_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, condu_ker%update(2,imate))
                 end if

              else if( words(1) == 'SPECI' ) then
                 !
                 ! Specific heat
                 !
                 sphea_ker % kfl_exist                = 1
                 if (exists('MIXTU')) then
                    sphea_ker % wlaws(imate)             = 'CPMIX'
                 else if (exists('CPGPT')) then
                    sphea_ker % wlaws(imate)             = 'CPGPT'
                 else if (exists('TABLE')) then
                    sphea_ker % wlaws(imate)             = 'CPTAB'
                 else if (exists('CMCNP')) then
                    sphea_ker % wlaws(imate)             = 'CMCNP'
                    !
                 else if (exists(ker_polynomial_name)) then
                    call ker_polynomial_readat( sphea_ker, imate )
                    !
                 else
                    sphea_ker % wlaws(imate)             = words(2)
                 endif
                 sphea_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 sphea_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT SPECIFIC HEAT')

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, sphea_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, sphea_ker%update(2,imate))
                 end if

              else if( words(1) == 'DUMMY' ) then
                 !
                 ! Dummy
                 !
                 dummy_ker % kfl_exist                = 1
                 dummy_ker % wlaws(imate)             = words(2)
                 dummy_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 dummy_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT DUMMY')

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, dummy_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, dummy_ker%update(2,imate))
                 end if

              else if( words(1) == 'TURBU' ) then
                 !
                 ! Turbulent viscosity
                 !
                 turmu_ker % kfl_exist                = 1
                 if (exists('SMAGO')) then
                    turmu_ker % wlaws(imate)          = 'SMAGO'
                 else if (exists('WALE ')) then
                    turmu_ker % wlaws(imate)          = 'WALE '
                 else if (exists('VRMAN')) then
                    turmu_ker % wlaws(imate)          = 'VRMAN'
                 else if (exists('SIGMA')) then
                    turmu_ker % wlaws(imate)          = 'SIGMA'
                 else if (exists('ILSA ')) then
                    turmu_ker % wlaws(imate)          = 'ILSA '
                 else if (exists('SSTKO')) then
                    turmu_ker % wlaws(imate)          = 'SSTKO'
                 else if (exists('STDKE')) then
                    turmu_ker % wlaws(imate)          = 'STDKE'
                 else if (exists('SPALA')) then
                    turmu_ker % wlaws(imate)          = 'SPALA'
                 else if (exists('MIXIN')) then
                    turmu_ker % wlaws(imate)          = 'MIXIN'
                 else if (exists('KOMEG')) then
                    turmu_ker % wlaws(imate)          = 'KOMEG'
                 else
                    turmu_ker % wlaws(imate)          = words(2)
                 end if

                 turmu_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 turmu_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT TURBULENT VISCOSITY')
                 if (exists('LOGVA')) kfl_logva=1

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')

                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, turmu_ker%update(1,imate))
                    if( icha2 /= '     ' )  then
                       call whn2up(icha2, turmu_ker%update(2,imate))
                    end if
                 end if

             else if( words(1) == 'ABSOR' ) then
                 !
                 ! Absorption
                 !
                 absor_ker % kfl_exist                = 1
                 absor_ker % wlaws(imate)             = words(2)
                 absor_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 absor_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT DUMMY')

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ') then
                    call whn2up(icha1, absor_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, absor_ker%update(2,imate))
                 end if

              else if( words(1) == 'SCATT' ) then
                 !
                 ! Scattering
                 !
                 scatt_ker % kfl_exist                = 1
                 scatt_ker % wlaws(imate)             = words(2)
                 scatt_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 scatt_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT DUMMY')

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, scatt_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, scatt_ker%update(2,imate))
                 end if

              else if( words(1) == 'MIXIN' ) then
                 !
                 ! MIXING
                 !
                 mixin_ker % kfl_exist                = 1
                 mixin_ker % wlaws(imate)             = words(2)
                 mixin_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 mixin_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT DUMMY')

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, mixin_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, mixin_ker%update(2,imate))
                 end if

              else if( words(1) == 'ANIPO' ) then
                 !
                 ! Anisotropic porosity
                 !
                 anipo_ker % kfl_exist                = 1
                 anipo_ker % wlaws(imate)             = words(2)
                 anipo_ker % rlaws(1:mlapa_ker,imate) = param(3:2+mlapa_ker)
                 anipo_ker % value_default(imate)     = getrea('DEFAU',-1.0_rp,'#DEFAULT DUMMY')
                 anipo_ker % kfl_type                 = MATRIX_PROPERTY

                 icha1 = getcha('UPDP1','     ','#When2update')
                 icha2 = getcha('UPDP2','NEVER','#When2update')
                 if( icha1 /= '     ' ) then
                    call whn2up(icha1, anipo_ker%update(1,imate))
                    if( icha2 /= '     ' )  call whn2up(icha2, anipo_ker%update(2,imate))
                 end if

              end if

              call ecoute('ker_readat')
           end do

!
! ALFALFALF
! DEBUG DEBUG DEBUG
! all this will be removed once the new windkessel is programmed
!
!        else if( words(1) == 'WINDK') then
!           !
!           !.md<1>WINDK
!           !.md<field>WINDK
!           !.md<com>(PRESSURE WINDKESSEL MODEL) Pressure model to account for arterial networks
!           !.md<com>outside of the scenario.
!           !
!           !
!           cpres_ker(1) = getrea('CCON1',0.0_rp,'#C constant for cavity 1')
!           rpres_ker(2) = getrea('RCON1',0.0_rp,'#R constant for cavity 1')
!           cpres_ker(1) = getrea('CCON2',0.0_rp,'#C constant for cavity 2')
!           rpres_ker(2) = getrea('RCON2',0.0_rp,'#R constant for cavity 2')

        else if( words(1) == 'DENSI' ) then
           !
           ! Density
           !
           !.md<1>DENSITY = rp
           !.md<field>DENSITY
           !.md<com>Medium density
           denme = getrea('DENSI',0.0_rp,'#MEDIUM DENSITY')

        else if( words(1) == 'VISCO' ) then
           !
           ! Viscosity
           !
           !.md<1> VISCOSITY = rp
           !.md<field>VISCOSITY
           !.md<com>Medium viscosity
           visme = getrea('VISCO',0.0_rp,'#MEDIUM VISCOSITY')

        else if( words(1) == 'GRAVI' ) then
           !
           !.md<1>GRAVITY: NORM= real1, GX= real2, GY= real2, GZ= real3                         $ Gravity acceleration
           !.md<field>GRAVITY
           !.md<com>Define the gravity vector
           !.md<com>g = real1*(real2,real3,real4)/|(real2,real3,real4)|
           !
           ! Gravity
           !
           grnor    = getrea('NORM ',0.0_rp,'#Gravity norm')
           gravi(1) = getrea('GX   ',0.0_rp,'#x-component of g')
           gravi(2) = getrea('GY   ',0.0_rp,'#y-component of g')
           gravi(3) = getrea('GZ   ',0.0_rp,'#z-component of g')
           call vecuni(3_ip,gravi,dummr)

        else if( words(1) == 'THICK') then
           !
           !.md<1>THICKNESS= real                                                               $ Thickness of free surface
           !.md<field>THICKNESS
           !.md<com>The thickness is the size over which the properties
           !.md<com>are smoothed out. It is usually taken as the size of three elements in the
           !.md<com>zone of interest.
           !
           ! Interface thickness
           !
           thicl = getrea('THICK',0.0_rp,'#Interface thickness')

        else if( words(1) == 'GASCO' ) then
           !
           !.md<1>GAS_CONSTANT= real                                                            $ Gas constant (example: 8.3144621 J/K)
           !.md<field>GAS_CONSTANT
           !.md<com>gas constant to be used in perfect gas law.
           !.md<com>For example, R=8.3144621 J/K.
           !
           ! Gas constant R
           !
           gasco = getrea('GASCO',8.3144621_rp,'#GAS CONSTANT') ! Default in J/K

        else if( words(1) == 'UREFE' ) then
           !
           ! Reference velocity
           !
           u_ref = getrea('UREFE',-1.0_rp,'#Reference velocity')

        else if( words(1) == 'HREFE' ) then
           !
           ! Reference height
           !
           h_ref = getrea('HREFE',-1.0_rp,'#Reference height')

        else if( words(1) == 'KREFE' ) then
           !
           ! Reference roughness
           !
           k_ref = getrea('KREFE',-1.0_rp,'#Reference roughness')

        else if( words(1) == 'USREF' ) then
           !
           ! Reference ustar
           !
           usref = getrea('USREF',-1.0_rp,'#Reference ustar')

        else if( words(1) == 'GEOST' ) then
           !
           ! Geostrophic wind (CFDWind2 model)
           !
           if( nnpar > 3 .or. nnpar < 2 ) call runend('KER_READAT: WRONG GEOSTROPHIC WIND DIMENSION')
           do ipara = 1,nnpar
              windg(ipara) = param(ipara)
           end do

        else if( words(1) == 'ROUGH' ) then
           !
           !.md<1>ROUGHNESS: CONSTANT, VALUE= real/FIELD=int                                    $ Roughness of wall-like surfaces
           !.md<field>ROUGHNESS
           !.md<com>Roughness of the surfaces. This affects the wall law and therefore
           !.md<com>only boundaries with a wall law condition. The rougness can be a constant value
           !.md<com>or can be given node-wise by a field to be defined in *.dom.dat file.
           !
           ! Roughness
           !
           if( words(2) == 'CONST' ) then
              kfl_rough = 0
              rough_dom = getrea('VALUE',-1.0_rp,'#Roughness value')
           else if( words(2) == 'FIELD' ) then
              kfl_rough = getint('FIELD',1_ip,'#Roughness field function')
           end if

        else if( words(1) == 'CANOP' ) then ! canopy height
           !
           !.md<1>CANOPY_HEIGHT: CONSTANT     /FIELD=int     $ Canopy height
           !.md<field>CANOPY_HEIGHT
           !.md<com>CANOPY_HEIGHTCanopy height is used by porous material of type canopy
           !.md<com>It can be a constant value or can be given node-wise by a field
           !.md<com>to be defined in *.dom.dat file for nodes of canopy material.
           !
           ! Canopy height
           !
           if( words(2) == 'CONST' ) then
              kfl_canhe = 0
           else if( words(2) == 'FIELD' ) then
              kfl_canhe = getint('FIELD',-1_ip,'#Canopy height function')
           end if

        else if( words(1) == 'HEIGH' ) then ! height over terrain
           !
           !.md<1>HEIGHT_OVER_TERRAIN: WALLD  /FIELD=int     $ Height of a node over the terrain
           !.md<field>HEIGHT_OVER_TERRAIN
           !.md<com> Height over the terrain is used by porous material of type canopy.
           !.md<com>It can be given node-wise by a field to be defined in *.dom.dat file for nodes of canopy material
           !.md<com>or it can be approximated by the wall distance. The wall distance approximation is a poor one
           !.md<com>and should be removed in teh future. It is maintained to make it compatible with te initial canopy implementation.
           !
           ! Canopy height
           !
           if( words(2) == 'WALLD' ) then
              kfl_heiov = 0
           else if( words(2) == 'FIELD' ) then
              kfl_heiov = getint('FIELD',-1_ip,'#Height over the terrain function')
           end if

        else if( words(1) == 'LAD  ' ) then
           !
           !.md<1>LAD: CONSTANT     /FIELD=int     $ Canopy height
           !.md<field>LEAF AREA DENSITY
           !.md<com>LAD is used by porous material of type canopy
           !.md<com>It can be a constant value or can be given node-wise by a field
           !.md<com>to be defined in *.dom.dat file for nodes of canopy material.
           !
           ! Leaf Area density
           !
           if( words(2) == 'CONST' ) then
              kfl_canla = 0
           else if( words(2) == 'FIELD' ) then
              kfl_canla = getint('FIELD',-1_ip,'#Canopy LAD function')
           end if
        else if( words(1) == 'WALLL' ) then
           !
           !.md<1>WALL_LAW:  ABL/REICHARDT, VARIABLE [MULTIPLIER= real]/WALL_DISTANCE= real      $ Law of the wall and distance to the wall
           !.md<field>WALL_LAW
           !.md<com>Definition of the law of the wall. Two informations are required. First the kind of
           !.md<com>wall law, either ABL (atmospheric boundary layer) or REICHARDT (viscous, buffer and log layer):
           !.md<com>
           !.md<com>      -  Reichardt's law:
           !.md<com>       \f$ u^+ = 1/kap ln(1+0.4y+) + 7.8 [ 1 - exp(-y+/11) - y+/11 exp(-0.33 y+) ] \f$.
           !.md<com>      -  ABL wall law: u^+ = 1 / 0.41 ln[ (y+k0) / k0 ]
           !.md<com>
           !.md<com>Then, the distance of the boundary nodes (computational wall) to the physical wall must be prescribed.
           !.md<com>It can be constant or variable.
           !.md<com>In the VARIABLE case, the distance to the wall is computed as HALF the size of the first element in the normal direction  multiplied
           !.md<com>by the constant MULTIPLIER= real. That is, if you want wall distance equals size of first element use MULTIPLIER=2. This also works for exchange location.
           !.md<com>By default, the wall distance of wall node is zero, meaning that the
           !.md<com>computational and physical walls coincide.
           !
           !
           ! Wall law
           !
           if(exists('VARIA')) then
              kfl_delta = 1   ! Variable wall distance
              delmu_dom = getrea('MULTI',1.0_rp,'Multiplier factor for variable wall distance')
           else if(exists('WALL2')) then
              kfl_delta = -1  ! convert wall-law boundaries to no slip boundaries, this is an alternative to using a negative wall distance
           else if(exists('CONST')) then
              kfl_delta = 0   ! Constant wall distance
           end if

           if(exists('ABL  ')) then
              kfl_ustar = 1
           else if(exists('REICH')) then
              kfl_ustar = 0
           else if(exists('ABL2 ')) then !ABL imposing zero gradient of k at the wall
              kfl_ustar = 2         !And imposing ustar from k and vel at nastin
              Cmu_st =getrea('CMU  ',0.0_rp, 'Cmu constant in k turbulent equation')
              if (Cmu_st < 1.0e-6_rp) call runend('ker_readat: Cmu needs to be positive when imposing ABL2 wall law')
           else if(exists('GRADP')) then
              kfl_ustar = 3
           end if

           if(exists('WALLD')) then
              delta_dom = getrea('WALLD',0.0_rp,'#Distance to the wall')
           end if

        else if( words(1) == 'WALLE' ) then
           kfl_waexl_ker = 1_ip
           if (exists('IMPLI'))  kfl_waexl_imp_ker =1_ip
           !
           ! For the moment I leave it like this
           ! If you want the exchange position to be calculated only for boundaries with codes 16, 18 & 20 you put
           !   WALL_EXCHANGE_LOCATION , BOUN1=16 , BOUN2=18, , BOUN3=20
           ! Perhaps a better option in order to avoid beeing limited to only 3 codes would be to add END_WALLE
           ! And then inside it BOUND  ... END_BOUND
           !
           if (exists('BOUN1'))  then
              kfl_boexc_ker = 0
              iauxi = getint('BOUN1',mcodb+1_ip,'#boundary code where to obtain exchange point')
              if(iauxi > mcodb) call runend('ker_readat:iauxi > mcodb - value should be given')
              kfl_boexc_ker(iauxi) = 1
           end if
           if (exists('BOUN2'))  then
              iauxi = getint('BOUN2',mcodb+1_ip,'#boundary code where to obtain exchange point')
              if(iauxi > mcodb) call runend('ker_readat:iauxi > mcodb - value should be given')
              kfl_boexc_ker(iauxi) = 1
           end if
           if (exists('BOUN3'))  then
              iauxi = getint('BOUN3',mcodb+1_ip,'#boundary code where to obtain exchange point')
              if(iauxi > mcodb) call runend('ker_readat:iauxi > mcodb - value should be given')
              kfl_boexc_ker(iauxi) = 1
           end if

        else if( words(1) == 'TWOLA' ) then
           !
           !.md<2>TWO_LAYER_MODEL: LES/RANS
           !.md<field>TWO_LAYER_MODEL
           !.md<com>Introduces a two layer wall model for LES, that uses an auxiliary RANS simulation
           !.md<com>to calculate the traction at the wall. This is achieved through coupling the two
           !.md<com>simulations, where the LES provides the velocity at the "inner" boundary of the RANS
           !.md<com>which in turn calculates the traction at the "outer" boundary (wall) and feeds it to
           !.md<com>to the LES. The two different options correspond to the two counterparts of the coupling.
           !.md<com>The LES velocity can either be instantaneous (not recommended) or averaged, for a period
           !.md<com>given by PERIOD = real (in the LES simulation). Code 22 needs to be used on the LES walls.
           !
           if( words(2) == 'LES  ') then
              kfl_twola_ker = 1
              if( words(3) == 'PERIO') then
                 tlape_ker = param(3)
              end if
           else if( words(2) == 'RANS ') then
              kfl_twola_ker = 2
           end if

!        else if( words(1) == 'NOSLW' ) then
!           !
!           !.md<2>NOSL_WALL_LAW
!           !.md<field>NOSL_WALL_LAW
!           !.md<com>   Applies wall law using no slip plus increased viscosity in teh first element
!           !.md<com>   Boundaries must be marked with 23
!           !
!           kfl_noslw_ker = 1

        else if( words(1) == 'RESTR' ) then

           krestr_2_codno(1:9) = nint(param(1:9))
           nrestr_2_codno = 9
           nrest: do ii=1,9
              if ( krestr_2_codno(ii) == 0_ip ) then
                 nrestr_2_codno = ii - 1
                 exit nrest
              end if
           end do nrest

        else if( words(1) == 'TIMEA' ) then
           !
           !.md<1>TIME_AVERAGING:       PERIOD=
           !
           kfl_wlaav_ker = 1_ip
           kfl_aveme_ker = 3           ! default ! Average - average - this is the correct option according to Yang et al.
           if( words(2) == 'PERIO' ) then
              tpeav_ker = param(2)
           end if
           if (exists('AVINS')) kfl_aveme_ker = 1
           if (exists('INSAV')) kfl_aveme_ker = 2
           if (exists('AVAV ')) kfl_aveme_ker = 3
        else if( words(1) == 'TENDE' ) then
           !
           !.md<1>TENDENCIES:      ON/OFF
           !
           if ( exists('ON   ').or.exists('YES  ')) kfl_tendencies_ker = .true.

           if (kfl_tendencies_ker) then ! read
              call read_tendencies() ! in mod_ker_tendencies
           end if

        else if( words(1) == 'COUPL' ) then
           !
           !.md<1>COUPLING
           !.md<2>char1 char2 int                                                             $ Module1, Module2, coupling number; Module1 <= Module2
           !.md<2>...
           !.md<1>END_COUPLING
           !.md<field>COUPLING
           !.md<com>This field contains the coupling information between modules.
           !.md<com>char1 and char2 are the modules and int is a number indicating
           !.md<com>the coupling procedure if different couplings exist.
           !.md<com>For example, "char1 char2 2" means that char2 modifies char1 module
           !.md<com>(char1 <= char2) using the coupling procedure number 2.
           !
           ! Coupling between modules
           !
           call ecoute('ker_readat')
           do while( words(1) /= 'ENDCO' )
              imodu = idmod(words(1))
              jmodu = idmod(words(2))
              if( jmodu == -1 .or. jmodu == -1 ) call runend('WRONG MODULES FOR COUPLING')
              kfl_coupl(imodu,jmodu) = max(1_ip,int(param(2),ip))
              !kfl_coupl(jmodu,imodu) = max(1_ip,int(param(2)))
              if( exists('FIELD') ) then
                 kfl_cowhe(imodu,jmodu) = getint('FIELD',1_ip,'#Field where to couple')
              else if( exists('ONFIE') ) then
                 kfl_cowhe(imodu,jmodu) = getint('ONFIE',1_ip,'#Field where to couple')
              end if
              call ecoute('ker_readat')
           end do

        end if
        call ecoute('ker_readat')
        !
        !.mdEND_PHYSICAL_PROBLEM
        !.md</code>
        !
     end do
     !--><group>
     !-->       <groupName>NUMERICAL_TREATMENT</groupName>
     !-->       <subGroup>
     !-->           <inputLine>
     !-->               <inputLineName>PRE_PENALIZATION</inputLineName>
     !-->               <inputLineHelp>Penalize pressure in time. To be used when converging to a steady state when the pressure presents an odd-even decoupling in time. This option helps converging by adding a damping effect.</inputLineHelp>
     !-->               <inputElement>
     !-->                   <inputElementType>edit</inputElementType>
     !-->                   <inputLineEditName>VALUE</inputLineEditName>
     !-->                   <inputLineEditValue>5</inputLineEditValue>
     !-->               </inputElement>
     !-->           </inputLine>
     !-->       </subGroup>
     !-->   </group>
     !
     !.md<sec>
     !.md# Numerical Treatment
     !.mdDefine the numerical treatment of kermod module. This field contains some options
     !.mdlike mesh multiplication or about the global arrays governed by PDE's needed by different modules
     !.md(like the generalized distance to the wall), etc. The difference of KERMOD with other MODULES
     !.mdis that when arrays are computed as a PDE, the boundary conditions are prescribed in the
     !.mdnumerical treatment field and not in a separate field.
     !.md<code>
     !.md<0>NUMERICAL_TREATMENT
     !
     ! Reach the section
     !
     do while( words(1) /= 'NUMER' )
        call ecoute('ker_readat')
     end do
     call ecoute('ker_readat')
     do while( words(1) /= 'ENDNU' )

        if( words(1) == 'FINIT' ) then
           !
           ! Finite volume data
           !
           if( option('FINIT') ) then
              kfl_lface     = 1   ! List of global faces
              kfl_fv_data   = 1   ! FV data structures
              kfl_elm_graph = 1   ! Element graph
           end if

        else if( words(1) == 'VECTO' ) then
           !
           !.md<1>VECTOR_SIZE
           !.md<field>VECTOR_SIZE
           !.md<com>Size of vectors for element assemblies.
           !
           if( words(2) == 'ADJUS' ) then
              kfl_vector_size = -1
           else
              kfl_vector_size = getint('VECTO',1_ip,'#vector size')
           end if
        else if( words(1) == 'TEMVE' ) then
           !To activate vectorization in temper
           if( words(2) == 'ON' ) then
              kfl_temper_vect = 1_ip
#ifndef VECTOR_SIZE        
              call messages_live('TEMPER VECTORIZATION: ON, BUT NO VECTOR_SIZE DEFINED') 
#else              
              call messages_live('TEMPER VECTORIZATION: ON') 
#endif
           end if
        else if( words(1) == 'ADAPT' ) then
           !
           !.md<1>ADAPTIVE_MESH_REFINEMENT
           !.md<field>ADAPTIVE_MESH_REFINEMENT
           !.md<com>Options concerning the mesh.
           !
           call AMR_readat()          

        else if( words(1) == 'MESH ' ) then
           !
           !.md<1>MESH
           !.md<field>MESH
           !.md<com>Options concerning the mesh.
           !
           call ecoute('ker_readat')
           do while(words(1) /= 'ENDME' )

              if ( words(1) == 'LOCAL' ) then
                 !
                 !.md<1>LOCAL_BASIS, NUMBER=int
                 !.md<2>LOCAL_BASIS= int1, BASIS: CARTESIAN/CYLINDRIC/SPHERICAL/CONSTANT PARAMETERS= real1,real2,...
                 !.md<2>LOCAL_BASIS= int2, BASIS: CARTESIAN/CYLINDRIC/SPHERICAL/CONSTANT PARAMETERS= real1,real2,...
                 !.md<2> ...
                 !.md<field>LOCAL_BASIS
                 !.md<com>Local coordinate systems for boundary conditions.
                 !
                 num_lobas = getint('NUMBE',1_ip,'#Number of coordinate systems')
                 call domain_memory_allocate('LOBAS',NUMBER1=num_lobas)
                 call ecoute('ker_readat')
                 do while( words(1) /= 'ENDLO' )
                    if( words(1) == 'LOCAL') then
                       icosi = getint('LOCAL',1_ip,'#Number of coordinate systems')
                       if ( words(2) == 'BASIS' ) then
                          nposi = 4_ip
                          wcsys = getcha('BASIS','     ','#Local coordinate system')
                          if (      wcsys == 'CARTE' ) then
                             !.md<com> - Set `BASIS: CARTESIAN` (default) to define a cartesian coordinate system by three points using
                             !.md<com>   `PARAMETERS=`<i>real1,real2,real3,real4,real5,real6,real7,real8,real9</i>. The first point <i>c(real1,real2,real3)</i>
                             !.md<com>   corresponds to the center of the coordinate system; The second point <i>a(real4,real5,real6)</i> has to form
                             !.md<com>   the primary axis P; The last point has to form the secondary axis S <i>b(rea7,rea8,rea9)</i>.
                             !.md<com>   The normal axis N is calculated internally doing the vectorial product of N = P x S.
                             lobas(icosi) % type = LOCAL_BASIS_CARTESIAN
                          else if ( wcsys == 'CYLIN' ) then
                             !.md<com> - Set `BASIS: CYLINDRICAL` to define a cylindrical coordinate system by three points using
                             !.md<com>   `PARAMETERS=`<i>real1,real2,real3,real4,real5,real6,real7,real8,real9</i>. The first point <i>c(real1,real2,real3)</i>
                             !.md<com>   corresponds to the center; The second point <i>a(real4,real5,real6)</i> has to form the radial axis R;
                             !.md<com>   The last point <i>b(real7,real8,real9)</i> has to form the axial axis A; The tangencial axis T is calculated
                             !.md<com>   internally doing the vectorial product of T = A x R.
                             lobas(icosi) % type = LOCAL_BASIS_CYLINDRICAL
                          else if ( wcsys == 'SPHER' ) then
                             lobas(icosi) % type = LOCAL_BASIS_SPHERICAL
                          else if ( wcsys == 'CONST' ) then
                             !.md<com> - Set `BASIS: CONSTANT` to use Alya coordinate system based on exterior normal. No points has to be defined.
                             lobas(icosi) % type = LOCAL_BASIS_CONSTANT
                          end if
                          lobas(icosi) % param(1:ndime*3) = param(nposi:ndime*3+nposi-1)
                       end if
                    end if
                    call ecoute('ker_readat')
                 end do
                 !
                 !.md<1>END_LOCAL_BASIS
                 !
              else if( words(1) == 'MULTI' ) then
                 !
                 !.md<2>MULTIPLICATION= int                                                        $ Number of mesh multiplication levels. 0 means do not multiply
                 !.md<field>MULTIPLICATION
                 !.md<com>Number of mesh multiplication levels. 0 means do not multiply. Postprocess can be carried
                 !.md<com>out on original mesh or on the last mesh of the multiplication process. It is strongly
                 !.md<com>advised to use extrapolated boundary conditions on boundaries so that the boundary conditions
                 !.md<com>on nodes are uniquely defined.
                 !
                 ! Mesh division
                 !
                 if( nnpar == 0 ) then
                    call ecoute('ker_readat')
                    do while(words(1) /= 'ENDMU' )
                       if( words(1) == 'LEVEL' ) then
                          !
                          !.md<2>LEVEL= int                                                        $ Level of multiplication.
                          !.md<field>LEVEL
                          !.md<com>Level of multiplication.
                          !
                          ndivi = getint('LEVEL',1_ip,'#Mesh divisions')
                       else if( words(1) == 'CURVA' ) then
                          !
                          !.md<2>CURVATURE= int                                                    $ Multiply elements locating new nodes on a geometry provided as an element field.
                          !.md<field>CURVATURE
                          !.md<com>1 means multiply with curvature
                          !.md<com>0 means standard straight-sided multiplication
                          !
                          multiply_with_curvature = getint('CURVA',1_ip,'#Divide with curvature')
                       else if( words(1) == 'PARAL' ) then
                          !
                          !.md<2>PARALLELIZATION: COORDINATES | GLOBAL                             $ Parallelization method for mesh multiplication.
                          !.md<field>PARALLELIZATION
                          !.md<com> Parallelization method to reconstruct the interfaces between domains. `GLOBAL` uses global element numbering (Default);
                          !.md<com> `COORDINATES` uses the coordinate method. When using interface elements or the element normal method multiplications
                          !.md<com> `COORDINATES` has to be activated.
                          !
                          if( words(2) == 'COORD' ) then
                             kfl_mmpar = 1
                          else if( words(2) == 'GLOBA' ) then
                             kfl_mmpar = 0
                          else
                             call runend('KER_READAT: MESH MULTIPLICATION UNKNOWN OPTION FOR PARALLELIZATION')
                          end if
                       else if( words(1) == 'ELEME' ) then
                          !
                          !.md<2>ELEMENT_NORMAL: ON/OFF                                            $ Multiply elements using element normal direction.
                          !.md<field>ELEMENT_NORMAL
                          !.md<com> Muliplication only for `QUA04` and `HEX08` elements based on the element normal direction.
                          !.md<com> Mesh must be oriented following a constant stacking direction.
                          !
                          if( option('ELEME') ) then
                             kfl_elndi = 1
                          else
                             kfl_elndi = 0
                          end if
                       end if
                       call ecoute('ker_readat')
                    end do
                    !
                    !.md<2>END_MULTIPLICATION
                    !
                 else
                    ndivi = getint('MULTI',1_ip,'#Mesh divisions')
                 end if

              else if( words(1)=='ARRAY') then
                 !
                 ! Element data base as an array
                 !
                 if( option('ARRAY') ) kfl_data_base_array = 1

              else if( words(1)=='SAVEE') then
                 !
                 !.md<1>SAVE_ELEMENT_DATA_BASE = Yes/No                           $ Save element data base (shape function derivatives, Hessian, etc.)
                 !.md<field>SAVE_ELEMENT_DATA_BASE
                 !.md<com>To avoid recomputing, each time they are needed, the shape functions derivatives
                 !.md<com>Hessian, etc., the element data base can be saved. Advantage: the code is faster.
                 !.md<com>Drawback: it requires more memory
                 !
                 if(      words(2) == 'PEREL' ) then
                    kfl_savda = 1
                 else if( words(2) == 'PERTY' ) then
                    kfl_savda = 3
                 else if(     option('SAVEE') ) then
                    kfl_savda = 1
                 end if

              else if( words(1) == 'DIVIS' ) then

                 ndivi = getint('DIVIS',1_ip,'#Mesh divisions')

              else if( words(1) == 'CURVA' ) then
                 !
                 !.md<2>CURVATURE= int                                                       $ Multiply elements locating new nodes on a geometry provided as an element field.
                 !.md<field>CURVATURE
                 !.md<com>1 means multiply with curvature
                 !.md<com>0 means standard straight-sided multiplication
                 !
                 ! Mesh division
                 !
                 multiply_with_curvature = getint('CURVA',1_ip,'#Divide with curvature')

              else if( words(1) == 'EDGEE' ) then
                 !
                 ! Edge elements
                 !
                 if( option('EDGEE') ) kfl_edge_elements = 1

              else if( words(1) == 'ROTAT' ) then
                 !
                 !.md<2>ROTATION: ANGLE=real, X/YZ                                           $ Rotation around an axe
                 !.md<field>ROTATION
                 !.md<com>Rotate the geometry around axe X, Y or Z.
                 !
                 ! Rotation
                 !
                 rotation_angle = getrea('ANGLE',0.0_rp,'#ANGLE OF ROTATION OF THE MESH')
                 if(      exists('X    ') ) then
                    kfl_rotation_axe = 1
                 else if( exists('Y    ') ) then
                    kfl_rotation_axe = 2
                 else if( exists('Z    ') ) then
                    kfl_rotation_axe = 3
                 else
                    call runend('KER_READAT: WRONG ROTATION AXE')
                 end if

              else if( words(1) == 'ELEME' ) then
                 !
                 !.md<2>ELEMENTAL_TO_CSR: ON/OFF                                             $ Elemental to CSR format saved in an array
                 !.md<field>ELEMENT_TO_CSR
                 !.md<com>Alya saves an array to go from elemental matrices to CSR format faster.
                 !
                 if( option('ELEME') ) kfl_element_to_csr = 1

              else if( words(1) == 'COO  ' ) then
                 !
                 !.md<2>COO: ON/OFF                                                          $ Compute COO format
                 !.md<field>COO
                 !.md<com>Compute COO format if it should be used in the solvers.
                 !
                 if( option('COO  ') ) kfl_coo = 1

              else if( words(1) == 'ELL  ' ) then
                 !
                 !.md<2>ELL: ON/OFF                                                          $ Compute ELL format
                 !.md<field>ELL
                 !.md<com>Compute ELL format if it should be used in the solvers.
                 !
                 if( option('ELL  ') ) kfl_ell = 1

              else if( words(1) == 'COOFO' ) then
                 if( option('COOFO') ) kfl_coo = 1

              else if( words(1) == 'FULLR' ) then
                 !
                 !.md<2>FULL_ROW: ON/OFF                                                     $ Compute full row graph
                 !.md<field>FULL_ROW
                 !.md<com>This option should be activiated if an algebraic solver uses
                 !.md<com>the full row format (/= partial row format).
                 !
                 if( option('FULLR') ) kfl_full_rows = 1

              else if( words(1) == 'EXTEN' ) then
                 !
                 !.md<2>EXTEND_GRAPH: OFF/ON
                 !.md<field>EXTEND_GRAPH
                 !.md<com>This option should be activiated if RAS preconditioning is used.
                 !.md<com>Add the connection between interface nodes in the
                 !.md<com>          original graph of the matrix and reallocate memory
                 !.md<com>          This is necessay to avoid the following possibility,
                 !.md<com>          where nodes 1 and 2 are not in the graph of subdomain 2.
                 !.md<com>
                 !.md<com>         o----o----o----o
                 !.md<com>         |    |    |    |  subdomain 1
                 !.md<com>         o----1----2----o
                 !.md<com>              |    |
                 !.md<com>              o----o
                 !.md<com>
                 !.md<com>         o----1    2----o
                 !.md<com>         |    |    |    |  subdomain 2
                 !.md<com>         o----o----o----o
                 !
                 if( option('EXTEN') ) kfl_graph = 1

              else if( words(1) == 'GRAPH' ) then
                 !
                 !.md<2>GRAPHS: ELEMENT, EXTENDED
                 !.md<field>GRAPHS
                 !.md<com>Combines all graphs options into a single option.
                 !
                 if( exists('ELEME') ) kfl_elm_graph = 1
                 if( exists('EXTEN') ) kfl_graph     = 1

              else if( words(1) == 'FACES' ) then
                 !.md<2>FACES_LIST: ON
                 if( option('FACES') ) kfl_lface = 1

              else if( words(1) == 'RENUM' .and. words(2) == 'NODES' ) then
                 !
                 ! Renumbering
                 !
                 if( words(3) == 'OFF  ' .or.  words(3) == 'NO   ' ) then

                    kfl_renumbering_npoin = 0

                 else if( words(3) == 'METIS' ) then

                    kfl_renumbering_npoin = 1

                 else if( words(3) == 'SFC  ' ) then

                    kfl_renumbering_npoin = 2
                    if( exists('NUMBE') ) &
                         nsfc_renumbering_npoin = getint('NUMBE',256_ip,'#NUMBER OF BINS IN EACH DIRECTION')
                    if( exists('SIZE ') ) &
                         nsfc_renumbering_npoin = getint('SIZE ',256_ip,'#NUMBER OF BINS IN EACH DIRECTION')

                 else if( words(3) == 'CUTHI' ) then

                    kfl_renumbering_npoin = 3

                 end if

              else if( words(1) == 'RENUM' .and. words(2) == 'ELEME' ) then
                 !
                 ! Renumbering
                 !
                 if( words(3) == 'OFF  ' .or.  words(3) == 'NO   ' ) then

                    kfl_renumbering_nelem = 0

                 else if( words(3) == 'NODES' ) then

                    kfl_renumbering_nelem = 1

                 else if( words(3) == 'EDGES' ) then

                    kfl_renumbering_nelem = 2

                 else if( words(3) == 'CLASS' ) then

                    kfl_renumbering_nelem = 3

                 end if

              else if( words(1) == 'BINEL' ) then
                 !
                 ! Element bin
                 !
                 element_bin_boxes(1)   = getint('NUMBE',3_ip,'#NUMBER OF BINS IN EACH DIRECTION')
                 element_bin_boxes(2:3) = element_bin_boxes(1)
                 kfl_element_bin        = 1

              else if( words(1) == 'SUPPO' ) then
                 !
                 !.md<2>SUPPORT_GEOMETRY, NODES=int1, BOUNDARIES=int2                              $ Number of mesh multiplication levels. 0 means do not multiply
                 !.md<3>  COORDINATES
                 !.md<4>  ...
                 !.md<3>  END_COORDINATES
                 !.md<3>  BOUNDARIES
                 !.md<4>  ...
                 !.md<3>  END_BOUNDARIES
                 !.md<2>END_SUPPORT_GEOMETRY
                 !.md<field>SUPPORT_GEOMETRY
                 !
                 kfl_suppo = 1
                 npoin_mm  = getint('NODES',0_ip,'#NUMBER OF NODES OF THE SUPPORT GEOMETRY')
                 nboun_mm  = getint('BOUND',0_ip,'#NUMBER OF NODES OF THE BOUNDARIES GEOMETRY')
                 ktype     = 0
                 if( npoin_mm + nboun_mm == 0 ) call runend('KER_READAT: WRONG DIMENSIONS OF THE SUPPORT GEOMETRY')
                 call ker_memory(5_ip)

                 call ecoute('ker_readat')
                 do while( words(1) /= 'ENDSU' )

                    if( words(1) == 'COORD' ) then
                       !
                       ! Coordinate COORD_MM(:,:)
                       !
                       call ecoute('ker_readat')
                       do while( words(1) /= 'ENDCO')
                          ipoin = int(param(1),ip)
                          if( ipoin < 0 .or. ipoin > npoin_mm ) then
                             call runend('KER_READAT: WRONG SUPPORT GEOMETRY COORDINATES')
                          end if
                          coord_mm(1:ndime,ipoin) = param(2:1+ndime)
                          call ecoute('ker_readat')
                       end do

                    else if( words(1) == 'TYPES' ) then
                       !
                       ! Connectivity LTYPB_MM(:,:)
                       !
                       call read_domain_arrays_types(0_ip,nboun_mm,ktype,ltypb_mm,lexis_mm)

                    else if( words(1) == 'BOUND' ) then
                       !
                       ! Connectivity LNODB_MM(:,:)
                       !
                       call ecoute('ker_readat')
                       do while( words(1) /= 'ENDBO')
                          iboun = int(param(1))
                          if( iboun < 0 .or. iboun > nboun_mm ) then
                             call runend('KER_READAT: WRONG SUPPORT BOUNDARY CONNECTIVITY')
                          end if
                          if( nnpar-1 /= mnodb_mm ) call runend('KER_READAT: WRONG SUPPORT BOUNDARY CONNECTIVITY')
                          lnodb_mm(1:mnodb_mm,iboun) = int(param(2:1+mnodb_mm),ip)
                          call ecoute('ker_readat')
                       end do

                    else if( words(1) == 'ALGEB' ) then
                       !
                       ! Algebraic solver
                       !
                       solve_sol => solve(3:)
                       call reasol(1_ip)

                    end if

                    call ecoute('ker_readat')
                 end do

              else if( words(1) == 'DEFOR' ) then
                 !
                 !.md<2>DEFORMATION
                 !.md<3>  METHOD: 1/2/3/4/5/6
                 !.md<3>  STEPS = 1
                 !.md<3>  ALGEBRAIC_SOLVER
                 !.md<4>  ...
                 !.md<3>  END_ALGEBRAIC_SOLVER
                 !.md<2>END_DEFORMATION
                 !.md<field>DEFORMATION
                 !
                 kfl_defor = 1
                 call ecoute('ker_readat')
                 do while( words(1) /= 'ENDDE' )

                    if( words(1) == 'METHO' ) then
                       !
                       ! Method
                       !
                       call runend('KER_READAT: NOT CODED')

                    else if( words(1) == 'STEPS' ) then
                       !
                       ! Steps
                       !
                       deformation_steps = getint('STEPS',1_ip,'#NUMBER OF DEFORMAITON STEPS')

                    else if( words(1) == 'CODES' .and. exists('NODES') ) then
                       !
                       !.md<2>CODES, NODES
                       !.md<2>  ...
                       !.md<2>END_CODES
                       !
                       tncod => tncod_ker(3:)
                       call reacod(1_ip)

                    else if( words(1) == 'ALGEB' ) then
                       !
                       ! Algebraic solver
                       !
                       solve_sol => solve(4:)
                       call reasol(1_ip)

                    end if

                    call ecoute('ker_readat')
                 end do

              else if( words(1) == 'CONSI' ) then
                 !
                 ! Consistent mass matrix
                 !
                 if( option('CONSI') ) then
                    kfl_conma = 1
                 end if

              else if( words(1) == 'WEIGH' ) then
                 !
                 ! Weighted consistent mass matrix
                 !
                 if( exists('CONSI') ) kfl_conma_weighted = 1
                 if( option('WEIGH') ) kfl_conma_weighted = 1

              else if( words(1) == 'INVER' ) then
                 !
                 ! APPROXIMATE INVERSE METHOD (NEW STUFF)
                 !
                 if( option('INVER') ) then
                    kfl_approx_inv_mass = 1
                 end if

              else if( words(1) == 'MASS ' ) then
                 if(      exists('APPRO') ) then
                    !
                    ! Approx inverse
                    !
                    kfl_approx_inv_mass = 1
                    kfl_dmass           = 1
                    kfl_conma           = 1

                 else if( exists('WEIGH') ) then
                    !
                    ! Consistent weighted mass
                    !
                    kfl_conma_weighted = 1

                 else if( exists('CONSI') ) then
                    !
                    ! Consistent mass
                    !
                    kfl_conma          = 1

                 else if( exists('DIAGO') ) then
                    !
                    ! Diagonally scaled mass
                    !
                    kfl_dmass          = 1
                 end if

              end if

              call ecoute('ker_readat')

           end do
           !
           !.md<1>END_MESH
           !
        else if( words(1) == 'GROUP' ) then
           !
           !.md<1>GROUPS
           !.md<2>  ALGEBRAIC_SOLVER
           !.md<3>  ...
           !.md<2>  END_ALGEBRAIC_SOLVER
           !.md<1>END_GROUPS
           !.md<field>GROUPS
           !
           call ecoute('ker_readat')
           do while( words(1) /= 'ENDGR' )

              if( words(1) == 'CODES' .and. exists('NODES') ) then
                 !
                 !.md<2>CODES, NODES
                 !.md<2>  ...
                 !.md<2>END_CODES
                 !
                 tncod => tncod_ker(6:)
                 call reacod(1_ip)

              else if( words(1) == 'ALGEB' ) then
                 !
                 ! Algebraic solver
                 !
                 solve_sol => solve(6:)
                 call reasol(1_ip)

              end if

              call ecoute('ker_readat')
           end do

        else if( words(1) == 'MATRI' ) then
           !
           ! Matrices to be computed by kernel
           !
           call ecoute('ker_readat')
           do while(words(1) /= 'ENDMA' )

              if( words(1) == 'GRADI' ) then
                 if( option('GRADI') ) kfl_matrix_grad = 1
              end if
              call ecoute('ker_readat')

           end do

        else if( words(1) == 'DIREC' ) then
           !
           ! Direct solver option
           !
           if(      words(2) == 'ALYA ' ) then
              kfl_direct_solver = SOL_DIRECT_SOLVER_ALYA
           else if( words(2) == 'SKYLI' ) then
              kfl_direct_solver = SOL_DIRECT_SOLVER_SKYLINE_ALYA
           else if( words(2) == 'PASTI' ) then
              kfl_direct_solver = SOL_DIRECT_SOLVER_PASTIX
           else if( words(2) == 'MUMPS' ) then
              kfl_direct_solver = SOL_DIRECT_SOLVER_MUMPS
           else if(words(2) ==  'WSMP ' ) then
              kfl_direct_solver = SOL_DIRECT_SOLVER_WSMP
           else if(words(2) ==  'PWSMP ' ) then
              kfl_direct_solver = SOL_DIRECT_SOLVER_PWSMP
           else if(words(2) ==  'GPUQR' ) then
              kfl_direct_solver = SOL_DIRECT_SOLVER_GPUQR
           end if

        else if( words(1) == 'GRADI' ) then
           !
           !.md<1>GRADIENT_PROJECTION: OPEN/CLOSE                                               $ Projection of gradients using open or close rule
           !.md<field>GRADIENT_PROJECTION
           !.md<com>Gradients are computed using an L2 projection, by solving the diagonal system
           !.md<com>M gi = \int du/dxi v dx. The option determines if the diagonal mass matrix M
           !.md<com>and the RHS are computed using a close rule or an open rule. In the last case,
           !.md<com>the mass matrix is lumped.
           !
           ! Gradient projection
           !
           if( words(2) == 'OPEN ' .or. words(2) == 'OPENE' ) then
              kfl_grpro = 0
           else
              kfl_grpro = 1
           end if

        else if( words(1) == 'DUALT' ) then
           !
           !.md<1>DUAL_TIME_STEP_PRECONDITIONER FACTOR/OFF 10
           !.md<field>DUAL_TIME_STEP_PRECONDITIONER
           !
           kfl_duatss = 1
           if( words(2) == 'FACTO') then
              fact_duatss = getint('FACTO',10_ip,'#Factor for dual time step precond')
           else
              kfl_duatss = 0
           end if
           print *,'kfl_duatss,fact_duatss',kfl_duatss,fact_duatss

        else if( words(1) == 'ELSES' ) then
           !
           !.md<1>ELSEST
           !.md<field>ELSEST
           !.md<com> ON/OFF
           !.md<com>Element search strategy. Alya may need to find the host element of some points,
           !.md<com>as this is the case of WITNESS_POINTS. The element search strategy can be tuned
           !.md<com>according to the mesh.
           !.md<com>Two strategies are available. The bin strategy is adapted to isotropic meshes.
           !.md<com>It consists in covering the computational volume by a Cartesian grid, and fill
           !.md<com>each box/cell with the elements which embedding box crosses the cell. The oct-tree strategy is more
           !.md<com>adapted to anisotropic meshes and creates a tree-like structure of boxes
           !.md<com>to accelerate the element search. For the bin strategy, the number of boxes
           !.md<com>in each direction should be prescribed in NUMBER_BOXES. For the oct-tree strategy,
           !.md<com>the user must give the maximum number of nodes (MAXIMUM_NUMBER) to fill the boxes. If the number of nodes
           !.md<com>exceeds this figure, then the box is further divided. If not, the box is filled with
           !.md<com>the elements  which embedding box crosses the box. The TOLERANCE option gives the tolerance to
           !.md<com>accept or reject a host element. A velue of 0.01 means that a host element is accepted
           !.md<com>if the node isoparametric coordinates in the element are within 1% of the size of the element (=1).
           !
           ! Elsest
           !
           if( option_not_off('ELSES') ) kfl_elses=1

           call ecoute('ker_readat')
           do while( words(1) /= 'ENDEL')

              if( words(1) == 'STRAT' ) then
                 !
                 !.md<2>STRATEGY: BIN/OCT_TREE                                                     $ Bin or quad/oct tree strategy
                 !
                 if(      words(2) == 'BIN  ' ) then
                    ielse(8) = 0
                 else if( words(2) == 'KDTRE' ) then
                    ielse(8) = 2
                 else
                    ielse(8) = 1
                 end if

              else if( words(1) == 'ELEME' ) then
                 !
                 ! Save element bounding boxes                                                         $ Save element bounding boxes
                 !
                 if( option('ELEME') ) ielse(16) = 1

              else if( words(1) == 'NUMBE' ) then
                 !
                 !.md<2>NUMBER_BOXES: int, int, int                                                $ Number of boxes in x,y,z directions for bin strategy
                 !
                 do idime = 1,3
                    ielse(idime) = max(int(param(idime),ip),1_ip)
                 end do

              else if( words(1) == 'MAXIM' ) then
                 !
                 !.md<2>MAXIMUM_NUMBER: int                                                        $ Max number of nodes in a box.
                 !
                 ielse(9) = getint('MAXIM',1_ip,'#MAXIMUM NUMBER OF NODES PER QUAD/OCT')

              else if( words(1) == 'OUTPU' ) then
                 !
                 !.md<2>OUTPUT: YES/NO                                                             $ Postprocess the bin or oct-tree (in GiD format)
                 !
                 if( words(2) == 'ON   ' .or. words(2) == 'YES  ' ) then
                    if(exists('STATI')) ielse( 7) = lun_elsta_dom
                    if(exists('MESH ')) ielse(12) = lun_elmsh_dom
                    if(exists('RESUL')) ielse(13) = lun_elres_dom
                 end if

              else if( words(1) == 'TOLER' ) then
                 !
                 !.md<2>TOLERANCE= real                                                            $ Tolerance to accept or reject a host element
                 !
                 relse(1) = getrea('TOLER',1.0e-6_rp,'#TOLERANCE')

              else if( words(1) == 'DATAF' ) then
                 if( words(2) == 'LINKE' ) then
                    ielse(4) = 1
                 else if( words(2) == 'TYPE ' ) then
                    ielse(4) = 0
                 end if

              else if( words(1) == 'MESHE' ) then
                 ielse(5) = getint('MESHE',1_ip,'#MAXIMUM NUMBER OF MESHES')
              end if

              call ecoute('ker_readat')
           end do
           !
           !.md<1>END_ELSEST
           !
        else if( words(1) == 'WALLD' ) then
           !
           !.md<1>WALL_DISTANCE
           !.md<field>WALL_DISTANCE
           !.md<com>This field contains the options (solver and boundary conditions) for the computation of the distance to the wall.
           !.md<com>The wall distance is added to the wall distance to the physical wall set in the PHYSICAL_PROBLEM field.
           !.md<com>By default, this distance is zero, meaning that the physical and computational walls coincide.
           !
           ! Example:
           !  WALL_DISTANCE
           !    ALGEBRAIC_SOLVER    DEFLATED_CG, ITERA= 5000,TOLER= 1e-9, ADAPTIVE, RATIO=1e-9
           !    PRECONDITIONING     DIAGONAL
           !    CODES, BOUNDARIES
           !      3 1
           !    END_CODES
           !  END_WALL_DISTANCE

           if(exists('SEARC')) then
              kfl_walld =  2
           elseif(exists('FILES')) then
              kfl_walld =  3
              kfl_walld_field(1) = -getint('FIELI',1_ip,'#Field Number for wall point index')
              kfl_walld_field(2) = -getint('FIELC',1_ip,'#Field Number for wall point coordinates')
           else
              kfl_walld =  1
           endif
           solve_sol => solve(2:)
           call ecoute('ker_readat')
           do while( words(1) /= 'ENDWA' )
              if( words(1) == 'ALGEB' ) then
                 !
                 !.md<2>ALGEBRAIC_SOLVER
                 !.md<2>  INCLUDE ./sample-solver.dat
                 !.md<2>END_ALGEBRAIC_SOLVER
                 !
                 call reasol(1_ip)
              else if( words(1) == 'PRECO' ) then
                 call reasol(2_ip)
              else if( words(1) == 'CODES' .and. exists('NODES') ) then
                 !
                 !.md<2>CODES, NODES
                 !.md<2>  INCLUDE ./sample-codes-nodes.dat
                 !.md<2>END_CODES
                 !
                 if(exists('GEOME')) then
                    tgcod => tgcod_ker(2:)
                    call reacod(4_ip)
                 else
                    tncod => tncod_ker(2:)
                    call reacod(1_ip)
                 end if
              else if( words(1) == 'CODES' .and. exists('BOUND') ) then
                 !
                 !.md<2>CODES, BOUNDARIES
                 !.md<2>  INCLUDE ./sample-codes-boundaries.dat
                 !.md<2>END_CODES
                 !
                 tbcod => tbcod_ker(2:)
                 call reacod(2_ip)
              end if
              call ecoute('ker_readat')
           end do
           !
           !.md<1>END_WALL_DISTANCE
           !
        else if( words(1) == 'WALLN' ) then
           !
           !.md<1>WALL_NORMAL
           !.md<field>WALL_NORMAL
           !
           !Example:
           !  WALL_NORMAL
           !    ALGEBRAIC_SOLVER    DEFLATED_CG, ITERA= 5000,TOLER= 1e-9, ADAPTIVE, RATIO=1e-9
           !    PRECONDITIONING     DIAGONAL
           !    CODES, BOUNDARIES
           !      3 1
           !    END_CODES
           !  END_WALL_NORMAL
           kfl_walln =  1
           solve_sol => solve(5:)
           call ecoute('ker_readat')
           do while( words(1) /= 'ENDWA' )
              if( words(1) == 'ALGEB' ) then
                 !
                 !.md<2>ALGEBRAIC_SOLVER
                 !.md<2>  INCLUDE ./sample-solver.dat
                 !.md<2>END_ALGEBRAIC_SOLVER
                 !
                 call reasol(1_ip)
              else if( words(1) == 'PRECO' ) then
                 call reasol(2_ip)
              else if( words(1) == 'CODES' .and. exists('NODES') ) then
                 !
                 !.md<2>CODES, NODES
                 !.md<2>  INCLUDE ./sample-codes-nodes.dat
                 !.md<2>END_CODES
                 !
                 tncod => tncod_ker(5:)
                 call reacod(1_ip)

              else if( words(1) == 'CODES' .and. exists('BOUND') ) then
                 !
                 !.md<2>CODES, BOUNDARIES
                 !.md<2>  INCLUDE ./sample-codes-boundaries.dat
                 !.md<2>END_CODES
                 !
                 tbcod => tbcod_ker(5:)
                 call reacod(2_ip)
              end if
              call ecoute('ker_readat')
           end do
           !
           !.md<1>END_WALL_NORMAL
           !
        else if( words(1) == 'CUTEL' ) then
           !
           !.md<1>CUT_ELEMENTS: ON/OFF                                                          $ If the mesh has cut elements
           !.md<field>CUT_ELEMENTS
           !.md<com>Cut elements are elements who requires a special integration rule
           !.md<com>when they are crossed by an interface. This option should be activated when
           !.md<com>simulating cracks with SOLIDZ or when using IMMBOU module to conserve
           !.md<com>immersed boundary volumes.
           !
           ! Cut elements
           !
           if( words(2) == 'ON   ' .or. words(2) == 'YES  ' ) kfl_cutel = 1

        else if( words(1) == 'HANGI' ) then
           !
           ! Hanging nodes
           !
           if( words(2) == 'ON   ' .or. words(2) == 'YES  ' ) kfl_hangi = 1

        else if( words(1) == 'HESSI' ) then
           !
           !.md<1>HESSIAN: ON/OFF                                                               $ If the mesh has cut elements
           !.md<field>HESSIAN
           !.md<com>Put OFF if Hessian of function forms should NOT be taken into account. The computational cost of the Hessian
           !.md<com>can be very high wrt to the gain of accuracy of certain algorithms.
           !
           if( words(2) == 'ON   ' .or. words(2) == 'YES  ' ) then
              kfl_lapla = 1
           else
              kfl_lapla = 0
           end if

        else if( words(1) == 'ROUGH' ) then
           !
           !.md<1>ROUGHNESS_EXTENSION
           !.md<field>ROUGHNESS_EXTENSION
           !.md<com>Extension of the roughness on the wall to the volume. This can be useful if the roughness
           !.md<com>of the nearest wall node is needed for an the interior node.
           !
           ! Roughness extension
           !
           kfl_extro = 1
           solve_sol => solve(1:)
           call ecoute('ker_readat')
           do while( words(1) /= 'ENDRO' )
              if( words(1) == 'ALGEB' ) then
                 !
                 !.md<2>ALGEBRAIC_SOLVER
                 !.md<2>  INCLUDE ./sample-solver.dat
                 !.md<2>END_ALGEBRAIC_SOLVER
                 !
                 call reasol(1_ip)
              else if( words(1) == 'PRECO' ) then
                 call reasol(2_ip)
              else if( words(1) == 'CODES' .and. exists('NODES') ) then
                 !
                 !.md<2>CODES, NODES
                 !.md<2>  INCLUDE ./sample-codes-nodes.dat
                 !.md<2>END_CODES
                 !
                 if(exists('GEOME')) then
                    tgcod => tgcod_ker(1:)
                    call reacod(4_ip)
                 else
                    tncod => tncod_ker(1:)
                    call reacod(1_ip)
                 end if
              else if( words(1) == 'CODES' .and. exists('BOUND') ) then
                 !
                 !.md<2>CODES, BOUNDARIES
                 !.md<2>  INCLUDE ./sample-codes-boundaries.dat
                 !.md<2>END_CODES
                 !
                 tbcod => tbcod_ker(1:)
                 call reacod(2_ip)
              end if
              call ecoute('ker_readat')
           end do
           !
           !.md<1>END_ROUGHNESS_EXTENSION
           !
        else if( words(1) == 'NOSLW' ) then
           !
           !.md<1>NO_SL_WALL_LAW
           !.md<field>NO_SL_WALL_LAW
           !.md<com>boundary code for no slip wall law - initially I was using kfl_fixbo_nsi == 23 but I need it to be in kernel
           !
           ! no slip wall law
           !
           kfl_noslw_ker = 1
           kfl_waexl_ker = 1      ! I also set wall exchange location because I want to feed frivel with the wall exchange velocity
           call ecoute('ker_readat')
           do while( words(1) /= 'ENDNO' )
              if( words(1) == 'CODES' .and. exists('BOUND') ) then
                 !
                 !.md<2>CODES, BOUNDARIES
                 !.md<2>  INCLUDE ./sample-codes-boundaries.dat
                 !.md<2>END_CODES
                 !
                 tbcod => tbcod_ker(7:)
                 call reacod(2_ip)
              end if
              call ecoute('ker_readat')
           end do
           !
           !.md<1>NO_SL_WALL_LAW
           !
        else if( words(1) == 'DISCR' ) then
           !
           !.md<1>DISCRETE_FUNCTIONS
           !
           call ker_discrete_function_read()
           !
           !.md<1>END_DISCRETE_FUNCTIONS
           !
           
        else if( words(1) == 'SPACE' ) then
           !
           !.md<1>SPACE_&_TIME_FUNCTIONS
           !
           call ecoute('ker_readat')
           do while( words(1) /= 'ENDSP' )
              if( words(1) == 'FUNCT' ) then

                 !-------------------------------------------------------------
                 !
                 ! Space/Time function definitions
                 !
                 !-------------------------------------------------------------
                 !
                 !.md<2>FUNCTION=char, DIMENSION=int
                 !.md<3>   f(x,t)
                 !.md<3>     ...
                 !.md<2>END_FUNCTION
                 !.md<field>SPACE_&_TIME_FUNCTIONS
                 !.md<com>Define space/time functions, depending on the coordinates x,y,z and time t.
                 !.md<com>SPACE_&_TIME_FUNCTIONS
                 !.md<com>   FUNCTION=RAMP [, FIELD]
                 !.md<com>      tanh(0.6*t) | NUMBER_FIELD= ifiel
                 !.md<com>   END_FUNCTION
                 !.md<com>END_SPACE_&_TIME_FUNCTIONS
                 !
                 number_space_time_function = number_space_time_function + 1
                 idime = 1
                 if( number_space_time_function > max_space_time_function ) call runend('KER_READAT: TOO MANY FUNCTIONS HAVE BEEN DEFINED')
                   if( exists('DIMEN') ) idime = getint('DIMEN',1_ip,'#Dimension of the function')
                   space_time_function(number_space_time_function) % ndime = idime
                   space_time_function(number_space_time_function) % name  = getcha('FUNCT','     ','#Space time function name')
                   igene = number_space_time_function
                   call ker_memory(7_ip)
                   do kdime = 1,idime
                      read(nunit,'(a)',err=3) space_time_function(number_space_time_function) % expression(kdime)
                      space_time_function(number_space_time_function) % expression(kdime) = adjustl(trim(space_time_function(number_space_time_function) % expression(kdime)))
                      !call ecoute('DONT')
                      !space_time_function(number_space_time_function) % expression(kdime) = adjustl(trim(ccard))
                      space_time_function(number_space_time_function) % nexpr = max( space_time_function(number_space_time_function) % nexpr , &
                           len(trim(space_time_function(number_space_time_function) % expression(kdime)),KIND=ip) )
                   end do
                   call ecoute('ker_readat')
                   if( words(1) /= 'ENDFU' ) call runend('KER_READAT: WRONG SPACE/TIME FUNCTION DIMENSION')
                end if
                call ecoute('ker_readat')
             end do
             !
             !.md<1>END_SPACE_&_TIME_FUNCTIONS
             !
          else if( words(1) == 'TIMEF' ) then

             !-------------------------------------------------------------
             !
             ! Time Function definitions
             !
             !-------------------------------------------------------------
             !
             !.md<1>TIME_FUNCTIONS
             !
             call ecoute('ker_readat')
             do while(words(1)/='ENDTI')
                !
                !.md<2>FUNCTION=char1, PARABOLIC | PERIODIC | DISCRETE, PARAMETERS= real1...real6      $ Time function definition
                !.md<3>...
                !.md<2>END_FUNCTION
                !.md<field>TIME_FUNCTIONS
                !.md<com>Define time function int1. Time functions are functions that multiply the values defined
                !.md<com>in the field "CODE, NODES", that is real1 real2 real3. Let t be the time
                !.md<com>and define te the effective time defined as te = min(max(t,real1),real2).
                !.md<com>The different options are:
                !.md<com>
                !.md<com>
                !.md<com>    PARABOLIC: f(t) = real3 * te^2 + real4 * te + real5.
                !.md<com>
                !.md<com>
                !.md<com>    PERIODIC: f(t) = real3 * cos(real4 * te + real5 ) + real6
                !.md<com>
                !.md<com>
                !.md<com>According to the definition of te, it means that:
                !.md<com>
                !.md<com>
                !.md<com>    f(t) = f(real1) for t < real1
                !.md<com>
                !.md<com>
                !.md<com>    f(t) = f(t) for real1 < t < real2
                !.md<com>
                !.md<com>
                !.md<com>    f(t) = f(real2) for t > real2
                !.md<com>
                !.md<com>
                !.md<com>For PARABOLIC and PERIODIC, no more informaiton is required and the function definition
                !.md<com>can be finalized by putting END_FUNCTION on the next line.
                !
                number_time_function = number_time_function + 1
                igene                = number_time_function
                ifunc                = number_time_function

                if( ifunc < 0 .or. ifunc > max_time_function ) then
                   call runend('ker_readat: WRONG FUNCTION NUMBER')
                end if

                if(      words(3) == 'PARAB' ) then
                   time_function(ifunc) % kfl_type = 1
                   time_function(ifunc) % npara    = 6

                else if( words(3) == 'PERIO' ) then
                   time_function(ifunc) % kfl_type = 2
                   time_function(ifunc) % npara    = 6

                else if( words(3) == 'CYCLI' ) then
                   time_function(ifunc) % kfl_type = 8
                   time_function(ifunc) % npara    = 2

                else if( words(3) == 'SNIF2' ) then
                   time_function(ifunc) % kfl_type = 9
                   time_function(ifunc) % npara    = 13

                else if( words(3) == 'DISCR' ) then
                   time_function(ifunc) % kfl_type = 3
                   time_function(ifunc) % npara    = 2 * getint('NUMBE',1_ip,'#Number of data')

                else if( words(3) == 'EAWAV' ) then
                   time_function(ifunc) % kfl_type = 10
                   time_function(ifunc) % npara    = 10

                else if( words(3) == 'KIAO ' ) then
                   time_function(ifunc) % kfl_type = 12
                   time_function(ifunc) % npara    = 0

                else

                   time_function(ifunc) % kfl_type = 0
                   number_time_function = number_time_function - 1

                end if

                if( time_function(ifunc) % kfl_type > 0 ) then
                   !
                   ! valid function
                   !
                   time_function(ifunc) % name = getcha('FUNCT','     ','#Time function name')
                   igene = ifunc
                   call ker_memory(8_ip)
                   if( time_function(ifunc) % kfl_type == 3 ) then
                      ifunp = 0
                      nfunp = time_function(ifunc) % npara / 2
                      if( nfunp < 1 ) call runend('KER_READAT: WRONG DISCRETE FUNCTION PARAMETER')
                      call ecoute('ker_readat')
                      do while( words(1) /= 'ENDFU' )
                         ifunp = ifunp + 1
                         if( ifunp > nfunp ) call runend('KER_READAT: WRONG DISCRETE FUNCTION DATA')
                         time_function(ifunc) % parameters( (ifunp-1)*2+1 ) = param(1)
                         time_function(ifunc) % parameters( (ifunp-1)*2+2 ) = param(2)
                         call ecoute('ker_readat')
                      end do
                      ! Order the function field
                      call ordena(nfunp,time_function(ifunc) % parameters)
                   else
                      dummi = time_function(ifunc) % npara
                      time_function(ifunc) % parameters(1:dummi) = param(4:dummi+3)
                   end if
                end if
                call ecoute('ker_readat')
                !
                !.md<1>END_FUNCTIONS
                !
             end do
          else if( words(1) == 'WINDK' ) then

             !-------------------------------------------------------------
             !
             ! Generic windkessel functions to use in Alya modules
             !
             !-------------------------------------------------------------
             !
             !.md<1>WINDKESSEL_FUNCTIONS
             !
             call ecoute('ker_readat')
             do while(words(1)/='ENDWI')
                !
                !.md<2>FUNCTION=char1  !!! ALF ALF ALF, PARABOLIC | PERIODIC | DISCRETE, PARAMETERS= real1...real6
                !.md<2>FUNCTION=char1  !!! ALF ALF ALF, TWO_ELEMENTS | THREE_ELEMENTS | FOUR_ELEMENTS, PARAMETERS= real1...real6
                !.md<3>...
                !.md<2>END_FUNCTION
                !.md<field>WINDKESSEL_FUNCTIONS
                !.md<com>Define Windkessel function char1. These should be included as boundary condition on boundaries (CODES, BOUNDARIES)
                !.md<com>The different options are:
                !.md<com>
                !.md<com>   TWO_ELEMENTS
                !.md<com>
                !.md<com>       q->
                !.md<com>         o----------+----------+
                !.md<com>                    |          |
                !.md<com>                    |          \
                !.md<com>                C1 ---      R1 /
                !.md<com>        p          ---         \
                !.md<com>                    |          /
                !.md<com>                    |          |
                !.md<com>         o----------+----------+
                !.md<com>
                !.md<com>
                !.md<com>   THREE_ELEMENTS
                !.md<com>
                !.md<com>       q->
                !.md<com>         o--/\/\/\-------+----------+
                !.md<com>             R2          |          |
                !.md<com>                         |          \
                !.md<com>                     C1 ---      R1 /
                !.md<com>        p               ---         \
                !.md<com>                         |          /
                !.md<com>                         |          |
                !.md<com>         o---------------+----------+
                !.md<com>
                !.md<com>   FOUR_ELEMENTS
                !.md<com>
                !.md<com>       q->        R2
                !.md<com>         o----+-/\/\/\-+-----+----------+
                !.md<com>              |        |     |          |
                !.md<com>              +--OOOO--+     |          \
                !.md<com>                 L1      C1 ---      R1 /
                !.md<com>        p                   ---         \
                !.md<com>                             |          /
                !.md<com>                             |          |
                !.md<com>         o-------------------+----------+
                !.md<com>
                !.md<com>
                !.md<com>
                !.md<com>
                !
                number_windk_systems = number_windk_systems + 1
                igene                = number_windk_systems
                ifunc                = number_windk_systems

                if( ifunc > max_windk_systems ) then
                   call runend('KER_READAT: MORE WINDKESSEL FUNCTIONS THAT ALLOWED BY MEMORY ALLOCATION')
                end if

                windk_systems(ifunc) % name = getcha('FUNCT','     ','#Time function name')

                if(      words(3) == 'ONEEL' ) then
                   windk_systems(ifunc) % wdks_model = 1
                   windk_systems(ifunc) % nparam     = 1
                   windk_systems(ifunc) % ndxs       = 1
                elseif(  words(3) == 'TWOEL' ) then
                   windk_systems(ifunc) % wdks_model = 2
                   windk_systems(ifunc) % nparam     = 2
                   windk_systems(ifunc) % ndxs       = 1
                else if( words(3) == 'THREE' ) then
                   windk_systems(ifunc) % wdks_model = 3
                   windk_systems(ifunc) % nparam     = 3
                   windk_systems(ifunc) % ndxs       = 1
                else if( words(3) == 'FOURE' ) then
                   windk_systems(ifunc) % wdks_model = 4
                   windk_systems(ifunc) % nparam     = 4
                   windk_systems(ifunc) % ndxs       = 2
                else
                  call runend('KER_READAT: TYPE OF WINDKESSEL FUNCTION NOT RECOGNISED')
                end if

                igene = ifunc
                call ker_memory(13_ip)
                dummi = windk_systems(ifunc) % nparam
                windk_systems(ifunc) % params(1:dummi) = param(3:dummi+2)

                call ecoute('ker_readat')

                !
                !.md<1>END_WINDKESSEL_FUNCTIONS
                !
             end do
          else if( words(1) == 'OPTIM' ) then
          !
          ! Optimization materials
          !                       kfl_opt_type => viscousity=1, conductivity=2, Arrhenius = 3
          !

             call ecoute('ker_readat')
             do while( words(1) /= 'ENDOP' )
                if( words(1) == 'OPTFU' ) then
                   kfl_cos_opt = 1_ip
                else if( words(1) == 'COSTT' ) then
                   kfl_cost_type = int(param(1))
                else if( words(1) == 'DVART' ) then
                   if (kfl_dvar_type /= 5) kfl_dvar_type = int(param(1))
                   !                 if (kfl_dvar_type == 5) kfl_ndvars_opt = ndime*npoin
                else if( words(1) == 'ADJOI' ) then
                   kfl_adj_prob = 1_ip
                else if( words(1) == 'NUMDE' ) then
                   kfl_ndvars_opt=int(param(1))
                endif
                call ecoute('ker_readat')
             end do
          else if (words(1) == 'REGUL') then
             kfl_regularization = .true.
             if( exists('SECON') ) kfl_second = .true.
             !
             !.md<1>REGULARIZATION SECOND EXP/EXP_LI/EXP_QU/EXP_GA/EXP_TR/EXP_CO/EXP_Q_C                                                               $ If the mesh has cut elements
             !

             if( exists('EXP  ') )  call set_regularization_type(reg_exp)
             if( exists('EXPLI') )  call set_regularization_type(reg_exp_linear)
             if( exists('EXPQU') )  call set_regularization_type(reg_exp_quadratic)
             if( exists('EXPGA') )  call set_regularization_type(reg_exp_garantzha)
             if( exists('EXPTR') )  call set_regularization_type(reg_exp_traslation)
             if( exists('EXPCO') )  call set_regularization_type(reg_exp_compression)
             if( exists('EXPQC') )  call set_regularization_type(reg_exp_quadratic_comp)
             if( exists('EXPLO') )  call set_regularization_type(reg_exp_log)
             if( exists('EXPCU') )  call set_regularization_type(reg_exp_cub)

             print *, 'REGULARIZATION, SECOND', kfl_regularization, kfl_second

          else if( words(1) == 'LOCAL' ) then
             !
             !.md<1>LOCALIZATION: UNILATERAL/BILATERAL     $ Contact type
             !.md<com>PDN contact algorithm
             !
             if( words(2) == 'UNILA') then
                kfl_conta = 1
             else if( words(2) == 'BILAT') then
                kfl_conta = 2
             end if

          else if( words(1) == 'LOOKU' ) then
             !
             !.md<1>LOOKUP_TABLES
             !
             !.md<field>LOOKUP_TABLES
             !.md<com>Define tables to be used for interpolating along arbitrary dimensions. 
             !.md<com>ORDER defines the order of the polynomial used for interpolation, 
             !.md<com>0 means snap on closest, 1 is linear, 2 is quadratic, etc.
             !.md<com>
             !.md<com>   E.g.:
             !.md<com>
             !.md<com>
             !.md<com>        TABLE: 1,   ORDER=1
             !.md<com>
             !.md<com>          INCLUDE table_from_FlameGen.dat
             !.md<com>
             !.md<com>        END_TABLE
             !.md<com>
             !.md<com>And define lookup frameworks 
             !.md<com>
             !.md<com>   E.g.:
             !.md<com>
             !.md<com>
             !.md<com>        FRAMEWORK: 1
             !.md<com>
             !.md<com>          MAIN_TABLE = 1
             !.md<com>
             !.md<com>          SCALING
             !.md<com>
             !.md<com>            ZMEAN  UNSCALED
             !.md<com>
             !.md<com>            ZVAR   VARIANCE [MEAN_INDEX=1]
             !.md<com>
             !.md<com>            CMEAN  TABLE 2
             !.md<com>
             !.md<com>            CVAR   OFF
             !.md<com>
             !.md<com>            IMEAN  TABLE 3
             !.md<com>
             !.md<com>          END_SCALING
             !.md<com>
             !.md<com>        END_TABLE
             !.md<com>
             !.md<com>   Or:
             !.md<com>
             !.md<com>        FRAMEWORK: 1
             !.md<com>
             !.md<com>          MAIN_TABLE = 1
             !.md<com>
             !.md<com>          SCALING
             !.md<com>
             !.md<com>            CMEAN  UNSCALED
             !.md<com>
             !.md<com>            CVAR   VARIANCE SQUARE
             !.md<com>
             !.md<com>            IMEAN  LIMIT  MIN=-123456.8  MAX=901.2
             !.md<com>
             !.md<com>          END_SCALING
             !.md<com>
             !.md<com>        END_TABLE
             !.md<com>
             !.md<com>
             !
             call ecoute('ker_readat')
             do while( words(1) /= 'ENDLO' )
                if( words(1) == 'TABLE' ) then
                   !.md<2>TABLE: int  [, ORDER=int]
                   !.md<3>  NVALU: int int int ...
                   !.md<3>  NVARI: int
                   !.md<3>  NCHUN: int 
                   !.md<3>  cordname1 SUBDI
                   !.md<3>   real1 real2 real3 ... real_NCHUN
                   !.md<3>   real  real ... real_NVALU1
                   !.md<3>  cordname2 SUBDI
                   !.md<3>   real1 real2 real3 ... real_NCHUN
                   !.md<3>   real  real ... real_NVALU2
                   !.md<3>  cordname3 SUBDI
                   !.md<3>  ...
                   !.md<3>  cordname1 cordname2 cordname3 ... varname1 varname2 varname3 ...
                   !.md<3>  real      real      real          real     real     real
                   !.md<3>  real      real      real          real     real     real
                   !.md<3>  real      real      real          real     real     real
                   !.md<3>  real      real      real          real     real     real
                   !.md<3>  real      real      real          real     real     real
                   i_lookup_tab = getint('TABLE',1_ip,'#Index of lookup table') 
                   if( i_lookup_tab > max_lookup_tab ) call runend('KER_READAT: lookup table index is over maximum')
                   number_lookup_tab = max(number_lookup_tab,i_lookup_tab)
                   order_loc = 1_ip
                   if (exists('ORDER')) then
                      order_loc = getint('ORDER',1_ip,'#Order of interpolation.')
                   endif
                   nullify(ptr_lookup_tab)
                   nullify(ptr_lookup_coords)
                   ptr_lookup_tab    => lookup_tab(i_lookup_tab)
                   ptr_lookup_coords => lookup_coords(:,i_lookup_tab)
                   call messages_live('LOOKUP TABLE '//trim(intost(i_lookup_tab))//', order of interpolation= '//trim(intost(order_loc)))
                   call tab_load_file(ptr_lookup_coords, ptr_lookup_tab, order=order_loc)
                   lookup_tab(i_lookup_tab) = ptr_lookup_tab
                   lookup_coords(:,i_lookup_tab) = ptr_lookup_coords
                   
                   do while(words(1)/='ENDTA')
                      call ecoute('ker_readat')
                   end do
                   !
                   !.md<2>END_TABLE
                   !.md<2>TABLE: int  [, ORDER=int]
                   !.md<2>...
                   !.md<2>END_TABLE
                   !
                elseif( words(1) == 'FRAME' ) then
                   !
                   !.md<2>FRAMEWORK: int
                   !
                   i_lookup_fw = getint('FRAME',1_ip,'#Index of lookup framework') 
                   if( i_lookup_fw > max_lookup_fw ) call runend('KER_READAT: lookup framework index is over maximum')
                   number_lookup_fw = max(number_lookup_fw,i_lookup_fw)
                   
                   call messages_live('LOOKUP FRAMEWORK '//trim(intost(i_lookup_fw)))
                   call messages_live('FRAMEWORK','START SECTION')
                   do while(words(1)/='ENDFR')
                      call ecoute('ker_readat')
                      if( words(1) == 'MAINT' ) then
                         !
                         !.md<3>MAIN_TABLE int                                                                $ index of main table
                         !
                         lookup_fw(i_lookup_fw) % kfl_tab_main =  getint('MAINT',1_ip,'#Index of main table')
                         lookup_fw(i_lookup_fw) % main_table   => lookup_tab(lookup_fw(i_lookup_fw) % kfl_tab_main)
                         call messages_live('Main table=                    ',INT_NUMBER = lookup_fw(i_lookup_fw) % kfl_tab_main,FMT='(a,1x,i5)')
                         
                         nullify(ptr_lookup_fw)
                         ptr_lookup_fw    => lookup_fw(i_lookup_fw)
                         call fw_allocate(1_ip, ptr_lookup_fw)
                         lookup_fw(i_lookup_fw) = ptr_lookup_fw
                      elseif( words(1) == 'SCALI' ) then
                         !
                         !.md<3>SCALING  
                         !
                         call ecoute('ker_readat')
                         do while(words(1)/='ENDSC')
                            !
                            !.md<4>cordname1: OFF | UNSCA | TABLE int | LIMIT MIN=real MAX=real | VARIA [MEANI=int] [SQUAR] | CONST=real 
                            !.md<4>cordname2: OFF | UNSCA | TABLE int | LIMIT MIN=real MAX=real | VARIA [MEANI=int] [SQUAR] | CONST=real
                            !.md<4>cordname3: OFF | UNSCA | TABLE int | LIMIT MIN=real MAX=real | VARIA [MEANI=int] [SQUAR] | CONST=real
                            !.md<4>...
                            !
                            
                            !
                            ! Find index of variable
                            !
                            ii = 0_ip
                            do jj = 1, lookup_fw(i_lookup_fw) % main_table % ndim
                                if (lookup_fw(i_lookup_fw) % main_table % coords(jj) % name == words(1)) ii = jj
                            enddo
                            
                            !
                            ! Handle error:
                            !
                            if (ii == 0) then
                                call runend('KER_READAT: '//words(1)//' is not a dimension of main table.')
                            endif

                            !
                            ! Assign scaling method
                            !
                            if( words(2) == 'OFF  ' ) then
                                lookup_fw(i_lookup_fw) % kfl_scale(ii) = -1
                                call messages_live('Scaling=                       ',CHAR_ARRAY =(/ words(1), 'OFF  ' /) ,FMT='(a,300(1x,a))')
                            elseif( words(2) == 'UNSCA' ) then
                                lookup_fw(i_lookup_fw) % kfl_scale(ii) = 0
                                call messages_live('Scaling=                       ',CHAR_ARRAY =(/ words(1), 'UNSCA' /) ,FMT='(a,300(1x,a))')
                            elseif( words(2) == 'TABLE' ) then
                                lookup_fw(i_lookup_fw) % kfl_scale(ii)          = 1
                                lookup_fw(i_lookup_fw) % scaling(ii) % kfl_tab  = getint('TABLE',1_ip,'#Index of scaling table')
                                icha1 = intost(lookup_fw(i_lookup_fw) % scaling(ii) % kfl_tab)
                                call messages_live('Scaling=                       ',CHAR_ARRAY =(/ words(1), 'TABLE', icha1 /) ,FMT='(a,300(1x,a))')
                            elseif( words(2) == 'LIMIT' ) then
                                lookup_fw(i_lookup_fw) % kfl_scale(ii)   = 2
                                lookup_fw(i_lookup_fw) % scaling(ii) % lims(1) = getrea('MIN  ',0.0_rp,'#Minimum for scaling')
                                lookup_fw(i_lookup_fw) % scaling(ii) % lims(2) = getrea('MAX  ',0.0_rp,'#Minimum for scaling')
                                icha1 = retost(lookup_fw(i_lookup_fw) % scaling(ii) % lims(1))
                                icha2 = retost(lookup_fw(i_lookup_fw) % scaling(ii) % lims(2))
                                call messages_live('Scaling=                       ',CHAR_ARRAY =(/ words(1), 'LIMIT', icha1, icha2 /) ,FMT='(a,300(1x,a))')
                            elseif( words(2) == 'CONST' ) then
                                lookup_fw(i_lookup_fw) % kfl_scale(ii)   = 3
                                lookup_fw(i_lookup_fw) % scaling(ii) % lims(1) = getrea('CONST',1.0_rp,'#Fixed scaled progress variable')
                                lookup_fw(i_lookup_fw) % scaling(ii) % lims(2) = lookup_fw(i_lookup_fw) % scaling(ii) % lims(1)
                                icha1 = retost(lookup_fw(i_lookup_fw) % scaling(ii) % lims(1))
                                call messages_live('Scaling=                       ',CHAR_ARRAY =(/ words(1), 'CONST', icha1 /) ,FMT='(a,300(1x,a))')
                            elseif( words(2) == 'VARIA' ) then
                                lookup_fw(i_lookup_fw) % kfl_scale(ii)   = getint('MEANI',ii-1_ip,'#Index of mean quantitiy corresponding to variance')
                                icha1 = intost(lookup_fw(i_lookup_fw) % kfl_scale(ii))
                                lookup_fw(i_lookup_fw) % kfl_scale(ii)   =  lookup_fw(i_lookup_fw) % kfl_scale(ii) +100 
                                icha2 = '     '
                                if (exists('SQUAR')) then
                                    lookup_fw(i_lookup_fw) % kfl_scale(ii)   = -1 * lookup_fw(i_lookup_fw) % kfl_scale(ii)
                                    icha2 = 'SQUAR'
                                endif
                                call messages_live('Scaling=                       ',CHAR_ARRAY =(/ words(1), 'VARIA', 'MEAN:', icha1, icha2  /) ,FMT='(a,300(1x,a))')
                            else
                                call runend('KER_READAT: unrecognised scaling option: '//words(2)//' of dimension: '//words(1))
                            endif

                            call ecoute('ker_readat')
                         end do
                         !
                         !.md<3>END_SCALING
                         !
                      endif
                   end do
                   call messages_live('FRAMEWORK','END SECTION')

                   !
                   ! Go through scaling dimensions and point to tables if necessary
                   !
                   do ii = 1, lookup_fw(i_lookup_fw) % main_table % ndim
                       select case (lookup_fw(i_lookup_fw) % kfl_scale(ii) )
                       case(1)
                           lookup_fw(i_lookup_fw) % scaling(ii) % tab => lookup_tab(lookup_fw(i_lookup_fw) % scaling(ii) % kfl_tab)
                       end select
                   enddo
                   nullify(ptr_lookup_fw)
                   ptr_lookup_fw    => lookup_fw(i_lookup_fw)
                   call fw_allocate(2_ip, ptr_lookup_fw)
                   lookup_fw(i_lookup_fw) = ptr_lookup_fw
                   !
                   !.md<2>END_FRAMEWORK
                   !.md<2>FRAMEWORK: int
                   !.md<2>...
                   !.md<2>END_FRAMEWORK
                   !
                endif
                call ecoute('ker_readat')
             enddo
             !
             !.md<1>END_LOOKUP_TABLES
             !
          end if
          call ecoute('ker_readat')
          !.mdEND_NUMERICAL_TREATMENT
          !.md</code>
       end do
       !--><group>
       !-->       <groupName>Output_postprocess</groupName>
       !-->       <groupType>oupos</groupType>
       !-->       <postProcess>
       !-->           <var>VELOCITY_MODULE</var>
       !-->           <var>VORTICITY</var>
       !-->           <var>KINETIC_ENERGY</var>
       !-->       </postProcess>
       !-->       <elementSet>
       !-->           <var>VELOCITY_MODULE</var>
       !-->           <var>VORTICITY</var>
       !-->           <var>KINETIC_ENERGY</var>
       !-->       </elementSet>
       !-->       <boundarySet>
       !-->           <var>MEAN_PRESSURE</var>
       !-->           <var>MASS</var>
       !-->           <var>FORCE</var>
       !-->           <var>TORQUE</var>
       !-->           <var>MEAN_YPLUS</var>
       !-->           <var>MEAN_VELOCITY</var>
       !-->           <var>WEAT_FORCE</var>
       !-->           <var>WEAT_SURFACE</var>
       !-->       </boundarySet>
       !-->       <nodeSet>
       !-->           <var>VELOX</var>
       !-->           <var>VELOY</var>
       !-->           <var>VELOZ</var>
       !-->           <var>PRESS</var>
       !-->           <var>YPLUS</var>
       !-->       </nodeSet>
       !-->       <witnessPoints>
       !-->           <var>VELOX</var>
       !-->           <var>VELOY</var>
       !-->           <var>VELOZ</var>
       !-->           <var>PRESS</var>
       !-->       </witnessPoints>
       !-->   </group>
       !.md<sec>
       !.md<0># Output and Post-process
       !.mdDefine the physical problem of kermod module. This field contains some variables shared by
       !.mdthe modules, like physical properties, distance to the wall, wall laws, etc.
       !.md<code>
       !.md<0>OUTPUT_&_POST_PROCESS
       !
       do while( words(1) /= 'OUTPU' )
          call ecoute('ker_readat')
       end do
       call ecoute('ker_readat')
       do while( words(1) /= 'ENDOU' )
          !
          !.md<1>POSTPROCESS char, STEPS= int                                                  $ Postprocess "char" variable each "int" time step
          !
          call posdef(2_ip,dummi)

          if( words(1) == 'VORTX' .or. words(1) == 'REDVE' ) then
             !
             kfl_vortx = 1
             if( words(2) == 'CELL ')  kfl_vortx = 1
             if( words(2) == 'FACE ')  kfl_vortx = 2
             if( words(3) == 'PARAL')  kfl_vortx_thres=2
             if( words(3) == 'THRES')  kfl_vortx_thres=1

          else if( words(1) == 'MESH ' ) then
             !
             ! Mesh output
             !
             if( option_not_off('MESH ') ) then
                kfl_oumes(1) = 1
                if( exists('GEOME') ) kfl_oumes(1) = 1
                if( exists('SETS ') ) kfl_oumes(2) = 1
                if( exists('BOUND') ) kfl_oumes(3) = 1
                if( exists('FIELD') ) kfl_oumes(4) = 1
                if( exists('ALL  ') ) kfl_oumes    = 1
             else
                kfl_oumes(1) = 0
             end if

          else if( words(1) == 'NOMES' ) then
             !
             ! Do not output mesh
             !
             kfl_oumes(1) = 0

          else if( words(1) == 'DETEC' ) then
             !
             ! Automatic detection
             !
             if( option('DETEC') ) then
                kfl_detection = 1
                call ecoute('ker_readat')
                do while( words(1) /= 'ENDDE' )
                   if( words(1) == 'LENGT' ) then
                      detection_length = getrea('LENGT',1.0_rp,'#Detection length')
                   else if( words(1) == 'VELOC' ) then
                      detection_velocity = getrea('VELOC',1.0_rp,'#Detection velocity')
                   end if
                   call ecoute('ker_readat')
                end do
             end if

          else if( words(1) == 'STL  ' ) then
             !
             ! Output boundary mesh STL format
             !
             if( words(2) == 'SUBDO' ) then
                kfl_oustl = 2
             else if( words(2) == 'BOUND' .or. words(1) == 'MASTE' ) then
                kfl_oustl = 1
             else if( words(2) == 'GID  ' ) then
                kfl_oustl = 3
             else
                if( words(2) /= 'OFF  ' .and.  words(2) /= 'NO   ' ) kfl_oustl = 1
             end if

          else if ( words(1) == 'THRES' ) then
             !
             ! Threshold
             !
             kfl_vortx_thres=1
             if( words(2) == 'VELOC' ) &
                  thr_veloc=getrea('VELOC',0.1_rp,'#velocity threshold')
             if( words(3) == 'LAMB2' ) &
                  thr_qvort=getrea('LAMB2',10.0_rp,'#lambda 2 threshold')

          else if( words(1) == 'ONLAS' ) then
             kfl_posdi = ndivi

          else if( words(1) == 'ONORI' ) then
             kfl_posdi = 0

          else if( words(1) == 'MESHQ' ) then
             !
             ! Mesh quality
             !
             if( words(2) == 'GAMMA' ) then
                kfl_quali = 1
             else if( words(2) == 'ASPEC' ) then
                kfl_quali = 2
             else if( words(2) == 'CONDI' ) then
                kfl_quali = 3
             else
                kfl_quali = 1
             end if

          else if( words(1) == 'ENSEM' ) then
             nsteps_ensemble = getint('ENSEM',0_ip,'#Number of steps for ensemble')

          else if( words(1) == 'STEPO' ) then
             npp_stepo = getint('STEPO',1_ip,'#Step over')
          else if( words(1) == 'STEPS' ) then
             npp_stepo = getint('STEPS',1_ip,'#Step over')

          else if( words(1) == 'FILTE' ) then
             !
             ! Filters
             !
             call fildef(2_ip)

          else if( words(1) == 'WITNE' .and. ( .not. exists('GEOME') ) .and. ( .not. exists('MESH ') ) .and. ( .not. exists('MESHE') ) ) then
             !
             !.md<1>WITNESS_POINTS, NUMBER=int, SCALE : XSCALE= real, YSCALE= real, ZSCALE= real
             !.md<1>...
             !.md<1>real1,real2,real3
             !.md<1>...
             !.md<1>END_WITNESS_POINTS
             !.md<field>WITNESS
             !.md<com>give the list of coordinates of witness points
             !
             ! COWIT: Witness points
             ! If number of witness points not given by user, use default value
             !
             if( words(2) == 'NUMBE' ) &
                  mwitn = getint('NUMBE',1_ip,'#NUMBER OF WITNESS POINTS')
             if( exists('SCALE') ) then
                scawi(1) = getrea('XSCAL',1.0_rp,'#x-factor')
                scawi(2) = getrea('YSCAL',1.0_rp,'#y-factor')
                scawi(3) = getrea('ZSCAL',1.0_rp,'#z-factor')
             end if
             !if( exists('METHO') ) then
             !    if( getcha('METHO','NULL ','#Witness poitn methodology') == 'COUPL')
             call ecoute('ker_readat')
             if( words(1) /= 'ENDWI' ) then
                call domain_memory_allocate('COWIT_ORIGI')
                do while( words(1) /= 'ENDWI' )
                   if( nwitn >= mwitn ) then
                      call runend('KER_READAT: TOO MANY WITNESS POINTS')
                   else
                      nwitn = nwitn + 1
                      do idime = 1,3
                         cowit_origi(idime,nwitn) = scawi(idime)*param(idime)
                      end do
                   end if
                   nwitn_all = nwitn
                   call ecoute('ker_readat')
                end do
             end if

          else if( words(1) == 'WITNE' .and. exists('GEOME') ) then
             !
             !.md<1>WITNESS, GEOMETRIES, NUMBER=int
             !.md<1>...
             !.md<1>real1,real2,real3
             !.md<1>...
             !.md<1>END_WITNESS
             !.md<field>WITNESS, GEOMETRIES
             !.md<com>give the list of witness geometries options
             !
             if( exists('NUMBE') ) &
                  nwitg = getint('NUMBE',1_ip,'#NUMBER OF WITNESS GEOMETRIES')
             call ecoute('ker_readat')
             call domain_memory_allocate('GEWIT')
             call witness_geometry_initialization()
             iwitg = 0
             do while( words(1) /= 'ENDWI' )
                if( iwitg >= nwitg ) then
                   call runend('KER_READAT: TOO MANY WITNESS GEOMETRIES')
                else
                   iwitg = iwitg + 1
                   call output_postprocess_read_geometry(gewit(iwitg) % kfl_geometry,gewit(iwitg) % param)                   
                end if
                call ecoute('ker_readat')
             end do

           else if( words(1) == 'WITNE' .and. ( exists('MESH ') .or. exists('MESHE') ) ) then
             !
             !.md<1>WITNESS, MESHES, NUMBER=int
             !.md<1>...
             !.md<1>SPHERE:      RADIUS=real1, X=real1,Y=real2,Z=real3, NAME=MESH1
             !.md<1>BOX:         XMIN=real1, XMAX=real2, YMIN=real3, YMAX=real4, ZMIN=real5, ZMAX=real6, NAME=MESH1
             !.md<1>BOUNDARY:    NAME=MESH1
             !.md<1>ELEMENT_SET: int1, NAME=MESH1
             !.md<1>MATERIAL:    int1, NAME=MESH1
             !.md<1>... 
             !.md<1>END_WITNESS
             !.md<field>WITNESS, MESHES
             !.md<com>give the list of witness geometries options to define witness meshes
             !
             if( exists('NUMBE') ) &
                  nwith = getint('NUMBE',1_ip,'#NUMBER OF WITNESS GEOMETRIES')
             call ecoute('ker_readat')
             call domain_memory_allocate('WITNESS_MESH')
             call witness_mesh_initialization()
             iwitg = 0
             do while( words(1) /= 'ENDWI' )
                if( iwitg >= nwith ) then
                   call runend('KER_READAT: TOO MANY WITNESS GEOMETRIES')
                else
                   iwitg = iwitg + 1
                   call output_postprocess_read_geometry(witness_mesh(iwitg) % geom % kfl_geometry,witness_mesh(iwitg) % geom % param)
                   if( exists('NAME ') ) then
                      witness_mesh(iwitg) % name = getcha('NAME ','NULL ','#Witness mesh name')
                   else
                      witness_mesh(iwitg) % name = trim(title)//'-'//'witness_mesh-'//trim(intost(iwitg))
                   end if
                end if
                call ecoute('ker_readat')
             end do

         else if( words(1) == 'PIXEL' ) then
             !
             ! Pixel output
             !
             kfl_pixel = 1
             call ecoute('ker_readat')
             do while( words(1) /= 'ENDPI' )
                if( words(1) == 'VARIA' ) then
                   if(      words(2) == 'VELOC' ) then
                      variable_pixel = ID_VELOC
                   else if( words(2) == 'PRESS' ) then
                      variable_pixel = ID_PRESS
                   else if( words(2) == 'TEMPE' ) then
                      variable_pixel = ID_TEMPE
                   else
                      call runend('KER_READAT: PIXEL VARIABLE NOT CODED, EDIT POSSPM.F90')
                   end if
                else if( words(1) == 'PLANE' ) then
                   if(      words(2) == 'X    ' ) then
                      plane_pixel = 1
                   else if( words(2) == 'Y    ' ) then
                      plane_pixel = 2
                   else if( words(2) == 'Z    ' ) then
                      plane_pixel = 3
                   else
                      call runend('KER_READAT: WRONG PIXEL PLANE')
                   end if
                else if( words(1) == 'COORD' ) then
                   coord_plane_pixel = getrea('COORD',0.0_rp,'#pixel plane coordinate')
                else if( words(1) == 'FREQU' ) then
                   kfl_pixel = getint('FREQU',1_ip,'#pixel output frequency')
                else if( words(1) == 'NUMBE' ) then
                   number_pixel(1) = int(param(1),ip)
                   number_pixel(2) = int(param(2),ip)
                   if( number_pixel(2) == 0 ) number_pixel(2) = number_pixel(1)
                end if
                call ecoute('ker_readat')
             end do

          else if( words(1) == 'VOXEL' ) then
             !
             ! Voxels
             !
             kfl_abovx = 1
             call ecoute('ker_readat')
             do while( words(1) /= 'ENDVO' )
                if( words(1) == 'BOUND' ) then
                   if( exists('AUTOM') ) then
                      kfl_abovx  = 1
                   else
                      kfl_abovx  = 2
                      bobvx(1,1) = param(1)
                      bobvx(1,2) = param(2)
                      bobvx(1,3) = param(3)
                      bobvx(2,1) = param(4)
                      bobvx(2,2) = param(5)
                      bobvx(2,3) = param(6)
                   end if
                else if( words(1) == 'RESOL' ) then
                   resvx(1) = int(param(1),ip)
                   resvx(2) = int(param(2),ip)
                   resvx(3) = int(param(3),ip)
                end if
                call ecoute('ker_readat')
             end do

          else if( words(1) == 'LIVEE' ) then
             !
             ! Give live info each [number] steps
             ! Live info is screen echo, convergence files, etc.
             ! Default is 1

             kfl_livee =  int(param(1))

          end if
          call ecoute('ker_readat')
       end do
       !
       !.md<0>END_OUTPUT_&_POST_PROCESS
       !.md</code>
       !
       if( kfl_timco == 2 ) then
          !        dtime = 1.0_rp
          timei = 0.0_rp
       end if
       !
       ! Types of Support surface have not been given
       !
       if( kfl_suppo == 1 ) then
          if( ktype == 0 ) then
             if( ndime == 2 ) then
                do iboun = 1,nboun_mm
                   ltypb_mm(iboun) = BAR02
                end do
             else
                do iboun = 1,nboun_mm
                   ltypb_mm(iboun) = TRI03
                end do
             end if
          end if
       end if
       if (kfl_waexl_ker/=0) then
          if (abs(delta_dom) < epsilon(1.0_rp) .and. kfl_delta /= 1 ) call runend('ker_readat:exchange location cannot be used with null wall distance ')
          ! When kfl_delta == 1 dexlo_ker will not be used - instead ywalb(iboun) will be used  - see ker_waexlo
          dexlo_ker =delta_dom ! set exchange locaton distance equal to wall distance
       end if
       !
       ! Kermod takes care of properties
       !
       kfl_prope = &
            &   densi_ker % kfl_exist  &
            & + visco_ker % kfl_exist  &
            & + poros_ker % kfl_exist  &
            & + condu_ker % kfl_exist  &
            & + sphea_ker % kfl_exist  &
            & + dummy_ker % kfl_exist  &
            & + turmu_ker % kfl_exist  &
            & + absor_ker % kfl_exist  &
            & + scatt_ker % kfl_exist  &
            & + mixin_ker % kfl_exist  &
            & + anipo_ker % kfl_exist
       !
       ! Velocity, temperature and concentration space time functions
       !
       if( kfl_vefun == 1001 ) then
          kfl_vefun = 1000 + space_time_function_number(wfname)
       end if
       if( kfl_tefun == 1001 ) then
          kfl_tefun = 1000 + space_time_function_number(wfnamt)
       end if
       if( kfl_cofun == 1001 ) then
          kfl_cofun = 1000 + space_time_function_number(wfnamc)
       end if
       if( kfl_difun == 1001 ) then
          kfl_difun = 1000 + space_time_function_number(wfnamd)
       end if

    end if
    !
    ! Load function number for subdomains
    !
    call ker_subdomain_function_name_to_number()

    return

3   call runend('KER_READAT: ERROR READING SPACE TIME FUNCTION')

  contains
    
  subroutine whn2up(ichari, update)
    !-------------------------------------------------------------------
    ! Routine that introduces the 'update' argument in prope_ker
    !
    !--------------------------------------------------------------------
    use def_master
    implicit none

    character(5), intent(in)   :: ichari
    integer(ip), intent(inout) :: update

    if (ichari =='ENDIN') then
       update = ITASK_ENDINN
    elseif (ichari =='ENDIT') then
       update = ITASK_ENDITE
    elseif (ichari =='DOITE') then
       update = ITASK_DOITER
    elseif (ichari =='ENDST') then
       update = ITASK_ENDSTE
    elseif (ichari =='ENDBL') then
       update = ITASK_CONBLK
    elseif (ichari =='NEVER') then
       update = 0
    else
       call runend('ker_readat: nonvalid update property in ker.dat: '//ichari)
    end if

  end subroutine whn2up
  
end subroutine ker_readat
