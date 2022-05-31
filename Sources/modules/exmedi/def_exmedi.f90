!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    def_exmedi.f90
!> @date    12/03/2013
!> @author  Mariano Vazquez
!> @brief   Def file
!> @details Def file
!> @}
!------------------------------------------------------------------------
module def_exmedi

!-----------------------------------------------------------------------
!    
! Heading for the module subroutines
!
!-----------------------------------------------------------------------
  use def_kintyp
  use def_master, only : kfl_exm_max_nmaterials
!
! General
!------------------------------------------------------------------------
! Parameters
!------------------------------------------------------------------------
  integer(ip), parameter ::&
      lun_vinte_exm = 512
  real(rp), parameter :: &
       zeexm = epsilon(0.0_rp)
  integer(ip), parameter ::&                     ! # postprocess variables
       nvarp_exm=60, nvars_exm=20, nvecp_exm=20, nspep_exm=35, nscap_exm=60
  integer(ip), parameter ::&                     ! # max sets
       msets_exm=10
  integer(ip), parameter  ::&
       EXM_CELL_STEADY_VOLTAGE=0, &                          !variables to decide on steady state in cell model
       EXM_CELL_STEADY_CALCIUM=1
  integer(ip), parameter  ::&
       EXM_CELLTYPE_MAXID = 3_ip, &                              !max cell type ID
       EXM_CELLTYPE_EPI=3_ip, &
       EXM_CELLTYPE_ENDO=1_ip, &
       EXM_CELLTYPE_MID=2_ip                                     !3-Epicardium, 1-endocardium, 2-midmyocardium
!------------------------------------------------------------------------
! Physical problem: read, modified or defined in exm_reaphy FIRST CALL
!------------------------------------------------------------------------
!--BEGIN REA GROUP
  integer(ip) ::&
       kfl_timei_exm,&             ! Existence of dT/dt
       kfl_gemod_exm,&             ! General model type of ionic 
       kfl_voini_exm(kfl_exm_max_nmaterials),&         ! Initial voltage flag (per material)
       kfl_fento_exm,&             ! Specific sub model of fenton's model
       kfl_appli_exm,&             ! Existence of applied currents
       kfl_appty_exm,&             ! Type of applied currents
       kfl_appva_exm,&             ! What is applied: current, density current, voltage...
       kfl_stead_exm,&             ! Steady-state has been reached 
       kfl_cemod_exm,&             ! Cell model (propagation model): monodomain o bidomain
       kfl_ptrig_exm,&             ! Electrical stimuli starting either with time or triggered by pressure
       kfl_hfmod_exm(kfl_exm_max_nmaterials),&         ! Runs the model on MODIFIED mode
       kfl_stree_exm,&             ! Construct a streeter fiber model
       onecl_exm(2),&              ! Number of beats to steady state and Cycle length
       kfl_drugs_exm,&             ! includes changes due to drug administration
       kfl_heter_exm,&             ! heterogeneous cell model
       kfl_nodif_exm,&             ! when set to 1, do not compute diffusion terms
       kfl_drugsmate_exm(kfl_exm_max_nmaterials),&     ! includes changes due to drug administration
       kfl_hetermate_exm(kfl_exm_max_nmaterials),&     ! heterogeneous cell model
       kfl_user_specified_celltypes_exm(kfl_exm_max_nmaterials,EXM_CELLTYPE_MAXID), &          ! Flag to know whether to run initial condition ODEs for different cell types. Each element 0|1, for ENDO, MID, EPI.  1 by default
       kfl_ignore_steadystate_celltypes_exm(kfl_exm_max_nmaterials,EXM_CELLTYPE_MAXID), &     ! Flag to know whether to die if different cell types did not reach steady state. Each element 1(ignore)|0(terminate), for ENDO, MID, EPI. 0 by deafult
       kfl_steadystate_variable(kfl_exm_max_nmaterials), &     ! Which variable will be used to test for the steady state. See EXM_CELL_STEADY_* variables
       moneclmate_exm(2,kfl_exm_max_nmaterials),&      ! Number of beats to steady state and Cycle length
       kfl_atbhe_exm,&             ! includes apex to base heterogeneity
       kfl_hfmodmate_exm(kfl_exm_max_nmaterials),&     ! Runs the model on Heart Failure mode
       kfl_psecg_exm,&             ! Flag for pseudo_ecg
       kfl_paced_exm,&             ! Flag for pacing or not the hiPSC-CMs
       kfl_fract_diffusion_exm(kfl_exm_max_nmaterials),& ! Fractional diffusion activated (per material)
       kfl_inaga_exm(kfl_exm_max_nmaterials),&           ! Flag that indicates if Sodium current from original O'Hara is exchanged for Passini INa dynamics
       kfl_save_convergence_cellmodel,& !for ohara model: 0 - do not save, 1 - save files with sequences of currents and stuff
       kfl_save_init_cellmodel      ! dump initial conditions to hardcode cell model 0 - do not, 1 - dump

!  integer(ip),  pointer ::&

       
  integer(ip) ::         &
       ndofn_exm,      &         ! # of d.o.f. of the problem
       ndof2_exm,      &         ! ndofn_*ndofn_
       nmate_exm,      &         ! # of materials
       nevat_exm,      &         ! Element matrix dim.=(ndime+1)*nnode
       nstim_exm,      &         ! Number of initial stimuli
       nstis_exm,      &         ! If needed, set number of initial stimuli
       modfi_exm,      &         ! Fiber model
       modor_exm(2),   &         ! Orhtotropic fiber model
       modce_exm,      &         ! Cell type
       modst_exm,      &         ! Stimuli field
       modab_exm,      &         ! Apex to base gradient field
       ncomp_exm,      &         ! Number of components of the fields
       nstrb_exm,      &         ! Number of sets to generate streeter fiber field 
       !ituss_exm ,     &         ! Cell type :  3-Epicardium, 1-endocardium, 2-midmyocardium
       strbo_exm(3),&              ! Boundaries for the streeter model
       nrootecg_exm,      &         ! Number of electrodes to be considered in the pseudo_ecg
       kcopeecg_exm                 ! Pseudo-ecg counter


  real(rp) ::&
       dtinv_exm    ,&              ! 1/dt
       xmccm_exm    ,&              ! Chi membrane capacitance
       aploo_exm(3)    ,&           ! Reset time (1), loop time (2) for applied stimuli loop and counter (3)
       xmccmmate_exm(kfl_exm_max_nmaterials)    ,&      ! Chi membrane capacitance
       voini_exm(kfl_exm_max_nmaterials),&              ! Initial voltage (per material)
       cleng_exm    ,&              ! Reference characteristic length       
       apval_exm(  15000),&         ! Applied current intensity  
       aplap_exm(  15000),&         ! Applied current time lapse 
       apcen_exm(3,15000),&         ! Applied current center 
       aprea_exm(  15000),&         ! Reach 
       scond_exm,&                  ! Surface conductivity coefficient
       react_exm,&                  ! Reaction term
       gdiff_exm(2,3,kfl_exm_max_nmaterials),&          ! Global diffusivities
       xmopa_exm(20,kfl_exm_max_nmaterials),&           ! Physical model parameters
       ttparmate_exm(3,13,kfl_exm_max_nmaterials),&     ! parameters to identify normal vs heart failure 1-cell simulations
       drugdmate_exm(24,kfl_exm_max_nmaterials),&       ! parameters to introduce a drug effect on cell models
       vminimate_exm(EXM_CELLTYPE_MAXID,kfl_exm_max_nmaterials), &       ! new voltage(celltype,imate)
       ttpar_exm(3,12),&            ! parameters to identify normal vs heart failure 1-cell simulations
       drugd_exm(12),&              ! parameters to introduce a drug effect on cell models
       vmini_exm(2,3), &            ! new voltage
       poref_fhn_exm(2),&           ! Potin and potfi, reference potentials for the FHN model 
       xthri_exm,&                  ! Lower threshold for fiber gradient
       xthrs_exm,&                  ! Upper threshold for fiber gradient
       vauin_exm(31,3,kfl_exm_max_nmaterials),&          ! initial values of  vauxi for 3D simulation 
       vcoin_exm(14,3,kfl_exm_max_nmaterials),&          ! initial vaules of vconc for 3D simulation
       vaulo_exm(31,3),&            ! Auxiliary variables for single cell solution
       vcolo_exm(14,3),&            ! Concentrations for single cell solution
       viclo_exm(27,3),&            ! Currents for single cell solution
       fiaxe_exm(3),&               ! Ventricular axis
       elmlo_exm(3),&               ! New local voltage
       stran_endo_exm,&             ! Angle for streeter endo
       stran_epi_exm,&              ! Angle for streeter epi
       volcai_exm,&                 ! Volume integral of calcium
       coordecg_exm(3,256),&        ! Coordinates of the electrodes for the pseudo-ecg
       pseudecg_exm(256),&          ! ECG values for the pseudo-ecg
       frequecg_exm,&               ! ECG postprocess frequency (time)
       fract_diff_coef_exm(kfl_exm_max_nmaterials), &           ! Fractional diffusion coefficient S
       kfl_steadystate_tolerance(kfl_exm_max_nmaterials)        ! Tolerance to use to determine the steady state in the cell model. Leave -1 by default, then subroutine will determine the correct one.

  integer(ip) ::&
       fract_diff_nintp_exm(kfl_exm_max_nmaterials), &          ! Fractional diffusion number of integration points
       ngrou_exm,&                  ! # of cell currents groups
       nauxi_exm,&                  ! # of auxiliary variables
       nconc_exm,&                  ! # of concentration unknowns
       ngate_exm,&                  ! # of activation gate fields
       nvint_exm,&                  ! Volume integral postprocess frequency (time steps)
       nicel_exm                    ! # of cell ionic currents
  
!------------------------------------------------------------------------
! Numerical problem: read, modified or defined in exm_reanut
!------------------------------------------------------------------------


  integer(ip) ::&
       kfl_genal_exm,&              ! General algorithm type ---DEPRECATED---
       kfl_timet_exm,&              ! Time treatment
       kfl_tiacc_exm,&              ! Temporal accuracy
       kfl_ticel_exm,&              ! Cell model time advance method
       kfl_tisch_exm,&              ! Temporal accuracy
       kfl_goite_exm,&              ! Keep iterating
       kfl_shock_exm,&              ! Shock capturing type 
       kfl_comat_exm,&              ! Compute amatr
       kfl_weigh_exm,&              ! Weighting of dT/dt
       kfl_adres_exm,&              ! Subiterations adaptive flag.
       kfl_normc_exm,&              ! Norm of convergence
       kfl_algso_exm,&              ! Type of algebraic solver
       kfl_repro_exm,&              ! Stabilization based on residual projection
       kfl_nolim_exm,&              ! Non-linear correction method
       kfl_nolum_exm,&              ! Non-linear terms lumped or not
       kfl_gcoup_exm,&              ! Geometric coupling with SOLIDZ
       miinn_exm,&                  ! Max # of iterations
       msste_exm,&                  ! Time substepping number
       last_iters_exm,&             ! Solver iterations for the current sub-iteration
       mnoli_exm,&                  ! Max # of iterations for the non-linear intracellular problem
       msoit_exm,&                  ! Max # of solver iterations
       nkryd_exm,&                  ! Krylov dimension
       memor_exm(2),&               ! Memory counter
       itera_exm,&                  ! Internal iteration counter
       nunkn_exm

  real(rp) ::&
       dtcri_exm    ,&              ! Critical time step
       shock_exm    ,&              ! Shock capturing parameter
       sstol_exm    ,&              ! Steady state tolerance
       cotol_exm    ,&              ! Convergence tolerance
       corat_exm    ,&              ! Convergence tolerance ratio
       dtext_exm    ,&              ! Externally fixed time step
       safet_exm    ,&              ! Safety factor for time step
       solco_exm    ,&              ! Solver tolerance
       weigh_exm,    &              ! Weight of dU/dt in the residual
       tnoli_exm    ,&              ! Tolerance for the non-linear intracellular problem
       resid_exm(10),&              ! Residual for outer iterations
       resou_exm(10) ,&              ! Reference residual for outer iterations 
       resin_exm(10) ,&              ! Reference residual for inner (sub)iterations 
       resin_first_exm(10) ,&        ! Reference residual for first inner (sub)iterations 
       err01_exm(2),&               ! L1 error T
       err02_exm(2),&               ! L2 error T
       err0i_exm(2),&               ! Linf error T
       err11_exm(2),&               ! L1 error grad(T)
       err12_exm(2),&               ! L2 error grad(T)
       err1i_exm(2),&               ! Linf error grad(T)
       staco_exm(3),&               ! Stability constants
       cpu_exmed(2),&               ! CPU for the EXM problem
       cpu_ass_sol_exm(4)           ! CPU time assembly and solver at each iteration

!------------------------------------------------------------------------
! Physical problem: read, modified or defined in exm_reaphy SECOND CALL
!------------------------------------------------------------------------

!
! Physical properties used in the model
!

  integer(ip),  pointer ::&
       idima_exm(:),      &              ! Diagonal indices for amatr
       kgrfi_exm(:)                      ! Large gradient fibers label vector
       
  real(rp),     pointer ::&
       cedif_exm(:,:,:)  ,&             ! Ex/Intracellular Diffusivity (point-wise)
       grafi_exm(:,:,:)  ,&             ! Fiber orientation gradient (a tensor)
       fiber_exm(:,:)    ,&             ! Fiber
       sheet_exm(:,:)    ,&             ! ortho fiber 1
       normal_exm(:,:)   ,&             ! ortho fiber 2
       celty_exm(:,:)    ,&             ! Cell types 
       atbhe_exm(:,:)    ,&             ! Apex to base Conductance gradient 
       vdiag_exm(:)      ,&             ! Diagonal preconditioner 
!+MRV       
       fibe2_exm(:,:)    

!------------------------------------------------------------------------
! Boundary conditions: read, modified or defined  in exm_reabcs
!------------------------------------------------------------------------
!
! Boundary conditions
! 
  integer(ip)           ::   kfl_exboc_exm    ! Boundary conditions explicitly given
  type(bc_nodes), pointer              :: &     
       tncod_exm(:)                          ! Node code type
  type(bc_bound), pointer              :: &     
       tbcod_exm(:)                          ! Boundary code type

!------------------------------------------------------------------------
! Output and Postprocess: read, modified or defined  in exm_reaous
!------------------------------------------------------------------------

  real(rp)                 :: &
       thiso_exm(3)                 ! Isochrones trigger, threshold, auto flag    
!--END REA GROUP
!------------------------------------------------------------------------
! Derived variables (slave-local in MPI)
!------------------------------------------------------------------------
  logical                  :: &
       iloner,weparal
!
! Physical properties used in the model
!
  real(rp),     pointer ::&
       amatr_auxi_exm(:)            ! Auxiliary amatr

  real(rp),     pointer ::&
       appfi_exm(:)                 ! Applied current field 
  
  real(rp)                ::&
       vamxm_exm(2,10)              ! Max-min values for the variables

  integer(ip),  pointer ::&
       lmate_exm(:),&               ! Materials (element-wise)
       kfl_fixno_exm(:,:),&         ! Nodal fixity 
       kfl_fixbo_exm(:)             ! Element boundary fixity
  real(rp),     pointer ::&
       bvess_exm(:,:)               ! Essential bc (or initial) values 
  real(rp),     pointer ::&
       bvnat_exm(:,:,:)             ! Natural bc values
  !
  ! Cell model types parameters and variables
  !
  real(rp), pointer ::&
       vauxi_exm(:,:,:),&           ! Auxiliary variables
       ticel_exm(:),&               ! Total cell ionic current
       jicel_exm(:),&               ! Jacobian of the total cell ionic current
       vicel_exm(:,:,:),&             ! Individual Cell ionic currents
       qneto_exm(:)

  real(rp),       pointer ::&
       refhn_exm(:,:)               ! Recuperation potential (only) for FHN

  integer(ip),    pointer ::&
       lapno_exm(:),&                 ! Model materials 2 / applied currents (node-wise)
       kwave_exm(:)                   ! Detect a new depolarization wave. 0 downstroke was detected, 1 upstroke was detected

  integer(ip),  pointer :: isoch_modified(:) ! TRUE when isochrone was detected. Used as a trigger to save isochrones

  integer(ip) ::&
       kfl_refre_exm(10)              ! Compute the reference residual for the exm_cvgunk's itask  

  real(rp) :: &
       sms_conversion_currents_exm

    real(rp) :: &                     ! some local cputime measurement
         cpold_exm

!
! Module types
! 
  type modcon_sol
     !
     ! This type defines each of the model constants, for each current
     !
     integer(ip)         :: iacti         ! activation state (0 , 1)
     integer(ip)         :: ituss        ! tissue (1 ENDO, 2 EPI , 3 M Cell)
     real(rp)            :: value(10)     ! parameter value (until 10 parameters)
     character(len=8)    :: iname         ! current name
     character(len=12)   :: vanam(10)     ! parameter name
     
  end type modcon_sol

  type(modcon_sol)       :: cupar_exm(20,20)


contains
!!!!!

  function heavis(a,b)
    use def_kintyp
    implicit none
    real(rp) :: a,b,heavis
    
    heavis = 1.0_rp
    if (a<b) heavis = 0.0_rp
    
    
  end function heavis
 

end module def_exmedi
