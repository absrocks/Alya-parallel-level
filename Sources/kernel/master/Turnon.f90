!-----------------------------------------------------------------------
!> @addtogroup Turnon
!> @{
!> @file    Turnon.f90
!> @author  Guillaume Houzeaux
!> @brief   Turnon the run
!> @details Initial operatons: read general data, mesh and modules data
!>          Construct the mesh dependent arrays.
!>
!> @}
!-----------------------------------------------------------------------
subroutine Turnon

  use def_parame
  use def_elmtyp
  use def_master
  use def_inpout
  use def_domain
  use def_kermod
  use mod_finite_volume,        only : finite_volume_arrays
  use mod_bourgogne_pinotnoir
  use mod_unity_tests
  use mod_auto_tuning,          only : auto_tuning_SpMV_OpenMP
  use mod_messages,             only : messages_report
  use mod_getopt,               only : kfl_check_data_file
  use mod_communications,       only : PAR_INITIALIZE_NON_BLOCKING_COMM, PAR_BARRIER
  use mod_par_affinity,         only : par_affinity
  use mod_messages,             only : messages_live
  use mod_outfor,               only : outfor
  use mod_messages,             only : livinf
  use mod_parall_openmp,        only : parall_openmp_adjacency_ompss_unity_test
  use mod_alya2dlb,             only : alya2dlb_initialization
  
  implicit none

  integer(8) :: count_rate8
  real(rp)   :: time1,time2

  !---------------------------------------------------------------------------------
  !
  ! Start run, creates code communicators in case of coupling, and read problem data
  !
  !---------------------------------------------------------------------------------

  call bourgogne(1_ip)
  !
  ! Compute time rate
  !
  call system_clock(count_rate=count_rate8)
  rate_time = 1.0_rp / max(real(count_rate8,rp),zeror)
  !
  ! Splits the MPI_COMM_WORLD for coupling with other codes. This defines the Alya world
  !
  call par_code_split_universe()
  !
  ! Splits the MPI_COMM_WORLD for coupling
  !
  call par_code_split_world()
  !
  ! Initialize DLB, just after initializing MPI
  !
  call alya2dlb_initialization()
  !
  ! Read problem data
  !
  call reapro()
  !
  ! Initialize time counters
  !
  call setgts(1_ip)
  !
  ! Initializations
  !
  call inidom()
  !
  ! Define element types
  !
  call elmtyp()                            ! Element lists: NNODE, LTOPO, LDIME, LLAPL...
  !
  ! Open module data files
  !
  call moddef(1_ip)
  !
  ! Exceptional things with modules (IMMBOU)
  !
  call excmod(2_ip)
  !
  ! Enable and initialize non-blocking communications
  !
  call PAR_INITIALIZE_NON_BLOCKING_COMM()
  !
  ! Affinity of processes
  !
  call par_affinity()

  !----------------------------------------------------------------------
  !
  ! Read domain data
  !
  !----------------------------------------------------------------------

  call messages_live('READ MESH ARRAYS','START SECTION')
  call readom()
  call messages_live('READ MESH ARRAYS','END SECTION')

  !----------------------------------------------------------------------
  !
  ! Domain partition... and compute some useful domain parameters
  !
  !----------------------------------------------------------------------

  call PAR_BARRIER()
  call pardom()

  !----------------------------------------------------------------------
  !
  ! Export mesh
  !
  !----------------------------------------------------------------------
  
  call mpio_export_domain()

  !----------------------------------------------------------------------
  !
  ! Read kermod data... AND
  ! do some calculations required by Kermod= NMATE, NSUBD
  !
  !----------------------------------------------------------------------

  call domvar(0_ip)
  call Kermod(-ITASK_TURNON)
  if( kfl_check_data_file == 1) call runend('O.K.!')

  !----------------------------------------------------------------------
  !
  ! Create domain data
  !
  !----------------------------------------------------------------------

  call messages_live('CONSTRUCT DOMAIN','START SECTION')
  call domain()
  call messages_live('CONSTRUCT DOMAIN','END SECTION')

  call cputim(time1)

  !----------------------------------------------------------------------
  !
  ! Others
  !
  !----------------------------------------------------------------------

  call Adapti(ITASK_TURNON)

  !----------------------------------------------------------------------
  !
  ! Read module data
  !
  !----------------------------------------------------------------------

  call livinf( 3_ip,'NULL',zero)
  call moduls(ITASK_TURNON)
  !
  ! Output optimization options
  !
  call outfor(-73_ip,0_ip,' ')
  !
  ! Parallel stuffs: send some data computed from modules data
  !
  call Parall(48_ip) 

  call livinf(52_ip,'NULL',zero)
  ! 
  ! Hybrid parallelization
  !
  call par_element_loop()
  call par_boundary_loop()
!!$  BLOCK
!!$    use def_domain, only : ompss_domains
!!$    use mod_parall, only : list_elements_par
!!$    use def_domain, only : ompss_boundaries
!!$    use mod_parall, only : list_boundaries_par
!!$    use mod_parall, only : num_pack_nboun_par
!!$    use mod_parall, only : num_pack_par
!!$   call parall_openmp_adjacency_ompss_unity_test(&
!!$         ompss_domains,list_elements_par,num_pack_par,npoin,lnnod,lnods)
!!$    call parall_openmp_adjacency_ompss_unity_test(&
!!$         ompss_boundaries,list_boundaries_par,num_pack_nboun_par,npoin,lnnob,lnodb)
!!$  END BLOCK
  call outfor(-72_ip)

  !----------------------------------------------------------------------
  !
  ! Required arrays depending on modules data
  !
  !----------------------------------------------------------------------
  !
  ! List of required arrays:
  !
  call reqarr()
  !
  ! PELPO_2, LELPO_2, PELEL_2, LELEL_2 extended Graphs
  !
  call lelpo2()
  !
  ! LFCNT, LNCNT: Identify faces of contact elements
  !
  call cntelm()
  !
  ! ZDOM_*, R_DOM_*, C_DOM:*. Schur type solvers and Aii preconditioners. * = Aii, Aib, Abi, Abb
  !
  call solpre()
  !
  ! R_SYM, C_SYM: Symmetric graph if necessary
  !
  call symgra()
  !
  ! LELBF, LELFA: Global face graph and element boundary face
  !
  call lgface()
  !
  ! Hanging nodes
  !
  call hangin()
  !
  ! Element bin for neighboring elements
  !
  call elebin()
  !
  ! Finite volume arrays
  !
  call finite_volume_arrays(meshe(ndivi))

  !----------------------------------------------------------------------
  !
  ! Others
  !
  !----------------------------------------------------------------------
  !
  ! Allocate memory for filters and check them
  !
  call fildef(-1_ip)
  !
  ! Allocate memory for all unknowns of the problem and coefficients
  !
  !call cputim(time1)
  call memunk(1_ip)
  call Optsol(ITASK_TURNON)
  !
  ! Check warnings and errors
  !
  call outerr(1_ip)
  !
  ! Read restart file
  !
  call restar(1_ip)
  !
  ! Close module data files and open occasional output files
  !
  call moddef(2_ip)
  !
  ! Compute groups: UNDER DEVELOPEMENT
  !
  !!!call ker_groups()
  !
  ! Initialization of some variables and solvers
  !
  call inivar(2_ip)
  !
  ! Close data file/open occasional files
  !
  call openfi(3_ip)
  !
  ! Support geometry
  !
  call ker_submsh()
  !
  ! Auto-tuning operations
  !
  call auto_tuning_SpMV_OpenMP()
  !
  ! Write information
  !
  call outcpu_operations(meshe(ndivi))

  call outinf()
  call outlat(1_ip)
  !
  ! Report
  !
  call messages_report()

  call cputim(time2)
  cpu_start(CPU_ADDTIONAL_ARRAYS) = cpu_start(CPU_ADDTIONAL_ARRAYS) + time2 - time1
  !
  !call initialize coprocessing
  !
#ifdef CATA
  call coprocessorinitializewithpython("coproc.py",9)
  call livinf(202_ip,' ',0_ip)
#endif
   !
   ! Transient fields
   !
  call calc_kx_tran_fiel()

end subroutine Turnon
