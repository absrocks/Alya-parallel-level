!------------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @name    ToolBox for solvers
!> @file    mod_solver.f90
!> @author  Guillaume Houzeaux
!> @date    13/06/2014
!> @brief   ToolBox for solvers
!> @details ToolBox for solvers. Matrix can be one block or multiblock
!>
!>          Multiblock matrices are organized this way:
!>
!>           ndof1   ndof2
!>          <-----><------->
!>          +-----+--------+
!>          |     |        |
!>          | 1,1 |  1,2   |
!>          +-----+--------+
!>          |     |        |
!>          | 2,1 |  2,2   |
!>          |     |        |
!>          +-----+--------+
!>
!>          REMARK: subroutine solver(...) is in /kernel/master
!>
!>          The diffrenet solvers implemented in Alya:
!>
!>          solver     format      ja               matrix
!>          ----------------------------------------------------
!>          Alya       BCSR        local            square
!>          Maphys     COO         local            square
!>          AGMG       CSR         local+ordering   rectangular
!>          PSBLAS     COO         local+info_glob  rectangular
!>          AMGX       BCSR        global           rectangular
!>
!>          WSMP       CSR         global           rectangular
!>          Pastix     CSC                          rectangular
!>
!>
!>
!>
!> @{
!------------------------------------------------------------------------

module mod_solver
  use def_kintyp,         only :  ip,rp,lg
  use def_kintyp,         only :  soltyp
  use def_master,         only :  INOTMASTER
  use def_master,         only :  IMASTER
  use def_master,         only :  IPARALL
  use def_master,         only :  ISEQUEN
  use def_master,         only :  NPOIN_TYPE
  use def_master,         only :  current_zone,lninv_loc
  use def_master,         only :  npoi1,npoi2,npoi3
  use def_master,         only :  kfl_async
  use def_domain,         only :  nperi
  use def_domain,         only :  lperi
  use def_domain,         only :  lmast
  use def_domain,         only :  nzdom
  use def_domain,         only :  r_dom
  use def_domain,         only :  c_dom
  use def_domain,         only :  npoin
  use def_domain,         only :  nzone,nzdom
  use def_domain,         only :  vmass,vmasc
  use def_solver,         only :  SOL_SOLVER_RICHARDSON
  use def_solver,         only :  memit
  use def_solver,         only :  SOL_BCSR_FORMAT,SOL_COO_FORMAT
  use def_solver,         only :  SOL_NO_SOLVER,SOL_NODES
  use def_solver,         only :  SOL_NODES
  use def_solver,         only :  SOL_ELEMENTS
  use def_solver,         only :  SOL_MATRIX_HAS_CHANGED
  use def_solver,         only :  SOL_CSR_FORMAT
  use def_solver,         only :  SOL_COO_FORMAT
  use def_solver,         only :  SOL_ELL_FORMAT
  use def_solver,         only :  SOL_EDGES
  use def_solver,         only :  SOL_OMP_OFF 
  use def_solver,         only :  SOL_OMP_STATIC 
  use def_solver,         only :  SOL_OMP_GUIDED 
  use def_solver,         only :  SOL_OMP_DYNAMIC
  use def_coupli,         only :  ncoup_implicit_d
  use def_coupli,         only :  lcoup_implicit_d
  use mod_memory,         only :  memory_alloca
  use mod_memory,         only :  memory_deallo
  use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications, only :  PAR_INTERFACE_OWN_NODE_EXCHANGE
  use mod_communications, only :  PAR_INTERFACE_EDGE_EXCHANGE
  use mod_communications, only :  PAR_GHOST_ELEMENT_EXCHANGE
  use mod_communications, only :  PAR_SUM
  use mod_communications, only :  PAR_MAX
  use mod_graphs,         only :  graphs_permut_metis_postordering
  use mod_graphs,         only :  graphs_copyij
  use mod_matrix,         only :  matrix_diagonal_CSR
  use mod_matrix,         only :  matrix_diagonal_COO
  use mod_matrix,         only :  matrix_diagonal_ELL
  use mod_matrix,         only :  matrix_CSR_SpMV
  use mod_matrix,         only :  matrix_COO_SpMV
  use mod_matrix,         only :  matrix_ELL_SpMV
  use mod_matrix,         only :  matrix_assemble_element_matrix_to_CSR
  use mod_matrix,         only :  matrix_assemble_element_matrix_to_COO
  use mod_matrix,         only :  matrix_assemble_element_matrix_to_ELL
  use mod_direct_solver,  only :  direct_solver_type_initialization  
  use mod_couplings,      only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,      only :  I_AM_IN_COUPLING
  implicit none

  private

  real(rp), parameter :: zeror =  epsilon(1.0_rp)
  real(rp), pointer   :: solver_invdiag(:)
  
  interface solver_parallel_vector_L2norm
     module procedure solver_parallel_vector_L2norm_1,        &
          &           solver_parallel_vector_L2norm_2
  end interface solver_parallel_vector_L2norm

  interface solver_parallel_scalar_product
     module procedure solver_parallel_double_scalar_product,  &
          &           solver_parallel_double_scalar_product_2,&
          &           solver_parallel_scalar_product_1,       &
          &           solver_parallel_scalar_product_2_3,     &
          &           solver_parallel_scalar_product_2_4,     &
          &           solver_parallel_scalar_product_t     
  end interface solver_parallel_scalar_product
  
  interface solver_initialize_matrix_and_rhs
     module procedure solver_initialize_matrix_and_rhs_S,&
          &           solver_initialize_matrix_and_rhs_1
  end interface solver_initialize_matrix_and_rhs
  
  public :: solver_assemble_element_matrix        ! Assemble element matrix
  public :: solver_preprocess
  public :: solver_postprocess
  public :: solver_solve_system
  public :: solver_solve
  public :: solver_periodicity
  public :: solver_periodicity_rhs
  public :: solver_initialize_matrix_and_rhs
  public :: solver_parallel_scalar_product
  public :: solver_parallel_double_scalar_product
  public :: solver_parallel_vector_L2norm
  public :: solver_smvp
  public :: solver_initialization
  public :: solver_parallel_residual_and_norm     ! Parallel residual and residual norm
  public :: solver_array_residual_norm
  public :: solver_krylov_subspace_save
  public :: solver_krylov_subspace_initial_guess
  public :: solver_diagonal
  public :: solver_parallel_SpMV                  ! Parallel matrix vector product
  public :: solver_lumped_mass_system
  public :: solver_allocate_system_memory         ! Allocate memory for a system A,b
  public :: solver_condition_number               ! Compute the conditioning of a matrix
  public :: solver_define_groups_from_field       ! Compute groups from a field
  public :: solver_invdiag
  public :: solver_errors                         ! Detect errors
  public :: solver_SpMV
  public :: solver_impose_dirichlet_condition     ! Impose Dirichlet condition on the unknown

contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2016
  !> @brief   Solver initialization
  !> @details Initialize solver type
  !>
  !----------------------------------------------------------------------

  subroutine solver_initialization(solve)
    type(soltyp), intent(inout)  :: solve
    integer(ip)                  :: ii
    !
    ! Input
    !
    solve % kfl_algso           = SOL_NO_SOLVER  ! No solver: do not modify (used to check errors)
    solve % ndofn               = 1              ! D.o.f
    solve % kfl_recov           = 0              ! Recover original residual after slave exchange
    solve % kfl_ortho           = 1              ! Orthogonalization GMRES: Gramm-Schmidt=0, modified=1
    solve % kfl_limit           = 0              ! Algebraic limiter
    solve % kfl_adres           = 0              ! Adaptive residual
    solve % nrhss               = 1              ! # of simultaneous RHS
    solve % miter               = 1              ! Max. number of iterations: should be 1, for Richardson solvers
    solve % nkryd               = 10             ! Krylov dimensions
    solve % kfl_coarse          = 0              ! No coarse solver
    solve % solco               = 1.0e-6_rp      ! Solver tolerance
    solve % adres               = 0.1_rp         ! Adaptive residual tolerance
    solve % solmi               = 1.0e-6_rp      ! Minimum solver tolerance
    solve % kfl_cvgso           = 0              ! Convergence flag
    solve % lun_cvgso           = 0              ! Convergence unit
    solve % kfl_solve           = 0              ! Output flag
    solve % lun_solve           = 0              ! Output unit
    solve % lun_exsol           = 0              ! External solver output unit
    solve % kfl_penal           = 0              ! Penalization: method
    solve % penal               = 1.0_rp         ! Penalization: parameter
    solve % kfl_defpr           = 2              ! Multigrid preconditioner: smoother preconditioner
    solve % itpre               = 1              ! Preconditioner iterations
    solve % kfl_preco           = 0              ! Preconditioner
    solve % kfl_leftr           = 0              ! Left preconditioner
    solve % kfl_marke           = 0              ! Output in market format
    solve % kfl_force           = 0              ! Do not force solver continuity
    solve % kfl_format          = SOL_CSR_FORMAT ! Matrix format
    solve % kfl_where           = SOL_NODES      ! Where are the unknowns
    solve % kfl_clean_precond   = 1              ! Preconditioner/coarse solver is not always recomputed at each solve
    solve % kfl_normalization   = 0              ! Residual normalization strategy
    solve % kfl_block_ras       = 0              ! Full RAS
    solve % num_multiple_solves = 1              ! Multiple solves with the same matrix
    solve % kfl_full_rows       = 0              ! By default, matrix are partiall assembled
    solve % kfl_save_krylov     = 0              ! Krylov subspace saving
    solve % kfl_kappa           = 0              ! Condition number should not be computed
    solve % normalization       = 1.0_rp         ! Residual normalization value
    solve % conf_file           = ''             ! Configuration file for solvers
    solve % bnorm_min           = 1.0e-12_rp     ! Minimum norm of RHS to consider null system
    solve % kfl_roe_correction  = 0              ! No round off error corrections 
    solve % threshold           = -1.0_rp        ! Threshold for direct solver
    solve % omp_schedule        = SOL_OMP_STATIC ! OMP schedule
    solve % omp_chunk_size      = 1000           ! OMP chunk size
    solve % omp_interface       = 0              ! Do not use openMP for interface > 0 = chunk
    solve % kfl_iffix           = 0              ! Fixity not imposed by solver
    solve % kfl_dirichlet       = 1              ! Way to impose Dirichlet boundary conditions 
    !
    ! Output
    !
    solve % iters       = 0                     ! Number of iterations
    solve % itsol(1)    =  huge(1_ip)           ! Max. # of solver iterations: do not modify
    solve % itsol(2)    = -huge(1_ip)           ! Min. # of solver iterations: do not modify
    solve % itsol(3)    = 0                     ! Ave. # of solver iterations: do not modify
    solve % nsolv       = 0                     ! # of solves
    solve % resin       = 1.0_rp                ! Initial preconditioned residual
    solve % resfi       = 1.0_rp                ! Final preconditioned residual
    solve % resi2       = 1.0_rp                ! Initial non-preconditioned residual = ||b-Ax||/||b||
    solve % resf2       = 1.0_rp                ! Final non-preconditioned residual = ||b-Ax||/||b||
    solve % cputi       = 0.0_rp                ! Total CPU time
    solve % num_spmv    = 0                     ! Number of SMVP
    solve % num_dot     = 0                     ! Number of dot products
    solve % cpu_spmv    = 0.0_rp                ! CPU time of SMVP
    solve % cpu_dot     = 0.0_rp                ! CPU dot product
    solve % cpu_schur   = 0.0_rp                ! Schur complement solver timing
    solve % xorth       = 0.0_rp                ! Orthogonalty check for CG solver
    solve % kfl_exres   = 0                     ! Do not output exact residual
    solve % xdiag       = 1.0_rp                ! Coefficient for diagonal solvers
    solve % num_solves  = 0                     ! Number of solves
    solve % bnorm       = 0.0_rp                ! Norm of RHS
    solve % kappa       = -999.0_rp             ! Condition number of the matrix
    solve % lambda_min  = -999.0_rp             ! Minimum eigenvalue
    solve % lambda_max  = -999.0_rp             ! maximum eigenvalue
    !
    ! Derived parameters
    !
    solve % wprob       = 'PROBLEM NAME'        ! Should be defined in module xxx_inivar
    solve % wsolv       = 'SOLVER NAME'         ! Defined in soldef
    solve % wprec       = 'NONE'                ! Defined in soldef
    solve % nzmat       = 1                     ! Matrix size
    solve % nzmat_own   = 1                     ! Matrix size for full row matrix
    solve % nzrhs       = 0                     ! RHS size
    solve % nequa       = 0                     ! # nodes
    solve % nequa_own   = 0                     ! # own nodes
    solve % nequa_halo  = 0                     ! # total nodes (including halos)    ! new mumps
    solve % nunkn       = 0                     ! # unknowns (nequa*ndofn)
    solve % ncols       = 0                     ! # columns
    solve % ndof2       = 1                     ! D.o.f^2
    solve % nequ1       = 0                     ! Interior unknwons
    solve % nequ2       = 0                     ! Start boundary unknowns
    solve % nequ3       = 0                     ! End boundary unknowns
    solve % kfl_symme   = 0                     ! Solver has an unsymmetric assembly
    solve % kfl_symeq   = 0                     ! Equation is non-symmetric
    solve % kfl_cmplx   = 0                     ! Complex solver
    solve % kfl_assem   = SOL_MATRIX_HAS_CHANGED ! Matrix not assembled
    solve % kfl_schur   = 0                     ! Default solver is not of Schur type
    solve % kfl_schum   = 0                     ! If matrices come from a Schur complement
    solve % kfl_version = 0                     ! Old version
    solve % heade       = 0                     ! Header has not been written
    solve % fil_cvgso   = ' '                   ! File name
    solve % fil_solve   = ' '                   ! File name
    solve % nzpre       = 0                     ! Preconditioner # non-null coefficients
    solve % kfl_update_precond  = 1             ! Preconditioner must be updated
    solve % kfl_mask    = 0                     ! If mask should be used for dot product
    nullify(solve % mask)                       ! Mask
    !
    ! Matrix graph
    !
    solve % nnz = 0                             ! Size of matrix (number of non-zeros)
    solve % nnz_ell = 0                         ! Size of ELL matrix (number of non-zeros)
    nullify(solve % ia)                         ! CSR format: IA
    nullify(solve % ja)                         ! CSR format: JA
    nullify(solve % ia_full)                    ! CSR format: IA for full row format
    nullify(solve % ia_full_ini)                ! CSR format: IA for full row format, only own nodes
    nullify(solve % ia_full_end)                ! CSR format: IA for full row format, only other nodes
    nullify(solve % ja_full)                    ! CSR format: JA for full row format
    nullify(solve % rows)                       ! COO format: rows
    nullify(solve % cols)                       ! COO format: columns
    nullify(solve % cols_ell)                   ! ELL format: columns
    !
    ! Block Gauss-Seidel
    !
    solve % kfl_blogs       = 0                 ! Monolithic solution / Block Gauss-Seidel
    solve % nblok           = 0                 ! Gauss-Seidel number of blocks
    solve % ndofn_per_block = 1                 ! dof per block
    nullify(solve % lperm_block)
    nullify(solve % linvp_block)
    !
    ! Fixity and boudnary conditions
    !
    solve % kfl_bvnat       = 0                 ! Natural b.c. not imposed in solver
    nullify(solve % kfl_fixno)                  ! Fixity array
    nullify(solve % bvess)                      ! Prescribed values array
    nullify(solve % bvnat)                      ! Prescribed force array
    !
    ! Deflated
    !
    solve % ngrou           = -1                ! Deflated CG: # groups (automatic # groups)
    solve % nskyl           = 0                 ! Deflated CG: Skyline A' matrix size
    solve % icoml           = 0                 ! Deflated CG: Communication group
    solve % ifbop           = 0                 ! Deflated CG: If boundary conditions are imposed on boundary nodes
    solve % kfl_gathe       = 0                 ! Deflated CG: all reduce
    solve % kfl_defas       = 0                 ! Deflated CG: skyline
    solve % kfl_defso       = 0                 ! Deflated CG: drect solver
    nullify(solve % limpo)                      ! List of imposed nodes (deflated CG)
    nullify(solve % lgrou)                      ! Group list for deflated CG
    nullify(solve % iskyl)                      ! Deflated CG: Skyline index
    nullify(solve % idiag)                      ! Deflated CG: Pointer to skyline diagonal
    nullify(solve % iagro)                      ! Deflated CG: Sparse graph
    nullify(solve % jagro)                      ! Deflated CG: Sparse graph
    nullify(solve % lcoun)                      ! Deflated CG: All gather strategy
    nullify(solve % lbig)                       ! Deflated CG: All gather strategy
    nullify(solve % displ)                      ! Deflated CG: All gather strategy
    nullify(solve % disp4)                      ! Deflated CG: All gather strategy
    nullify(solve % lcou4)                      ! Deflated CG: All gather strategy
    nullify(solve % xsmall)                     ! Deflated CG: All gather strategy
    nullify(solve % xbig)                       ! Deflated CG: All gather strategy
    !
    ! Renumbered Gauss-Seidel
    !
    solve % kfl_renumbered_gs = 0               ! Renumbered Gauss-Seidel permutation
    nullify(solve % permr_gs)                   ! Renumbered Gauss-Seidel permutation
    nullify(solve % invpr_gs)                   ! Renumbered Gauss-Seidel inverse permutation
    nullify(solve % lgrou_gs)                   ! Renumbered Gauss-Seidel groups
    nullify(solve % vecto_gs)                   ! Renumbered Gauss-Seidel vetor
    nullify(solve % adiag1)                     ! Renumbered Bidiagonal storing subdiagonal terms
    nullify(solve % idiag1)                     ! Renumbered Bidiagonal storing the position of the column according to the subdiagona
    !
    ! Linelet
    !
    solve % kfl_factl = 0                       ! Linelet: not factorized
    solve % nline     = 0                       ! Linelet:
    solve % nlin1     = 0                       ! Linelet:
    solve % npoin     = 0                       ! Linelet: # of nodes in linelets
    solve % kfl_linty = 0                       ! Linelet: linelet type
    solve % toler     = 10.0_rp                 ! Linelet: tolerance aspect ratio. Small -> lots of linelet, large -> few
    nullify(solve % lpntr)                      ! Linelet: tridiag. to original sym. matrix
    nullify(solve % lrenu)                      ! Linelet: point position in tridiag.
    nullify(solve % lrenup)                     ! Linelet: inverse permutation
    nullify(solve % limli)                      ! Linelet: list of imposed node
    nullify(solve % lline)                      ! Linelet: Pointer for each linelet
    nullify(solve % trima)                      ! Linelet: Tridiagonal matrix
    !
    ! Schur solver
    !
    solve % kfl_scpre = 0                       ! Schur preconditioner solver: skyline
    solve % kfl_scaii = 1                       ! Schur matrix solver: sparse
    solve % poaii     = 0                       ! Pointer to Aii
    solve % poaib     = 0                       ! Pointer to Aib
    solve % poabi     = 0                       ! Pointer to Abi
    solve % poabb     = 0                       ! Pointer to Abb
    !
    ! Block structure
    !
    solve % num_blocks       = 1
    solve % block_num        = 1
    solve % block_dimensions = 1
    solve % block_pointers   = 1
    solve % kfl_react        = 0
    nullify( solve % reaction       )
    nullify( solve % lpoin_block    )
    nullify( solve % lpoin_reaction )
    nullify( solve % bvess          )
    nullify( solve % kfl_fixno      )
    nullify( solve % bvnat          )
    do ii = 1,size(solve % block_array,KIND=ip)
       nullify( solve % block_array(ii) % bvess     )
       nullify( solve % block_array(ii) % kfl_fixno )
       nullify( solve % block_array(ii) % bvnat     )
    end do
    !
    ! Main matrix comes from a Schur complement
    ! A = A1 + A2*A3^-1*A4
    !
    solve % ndofn_A3 = 1             ! Dimension of A3
    nullify(solve % A1)
    nullify(solve % A2)
    nullify(solve % A3)
    nullify(solve % A4)
    nullify(solve % invA3)
    !
    ! Krylov subspace savings
    !
    nullify(solve % krylov_dir)
    nullify(solve % krylov_dot_products)
    solve % krylov_size = 0
    !
    ! Direct solvers
    !
    call direct_solver_type_initialization(solve % direct_solver)           ! Direct solver
    nullify(solve % direct_solver_RAS)                                      ! Direct solver for RAS
    call direct_solver_type_initialization(solve % direct_solver_coarse)    ! Direct solver for coarse
    call direct_solver_type_initialization(solve % direct_solver_block_LU)  ! Direct solver for block LU
    call direct_solver_type_initialization(solve % direct_solver_AMG)       ! Direct solver for AMG
    call direct_solver_type_initialization(solve % direct_solver_Deflation) ! Direct solver for deflation
 
  end subroutine solver_initialization

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/05/2017
  !> @brief   Save a Krylov subspace
  !> @details Save a Krylov subspace
  !>
  !----------------------------------------------------------------------

  subroutine solver_krylov_subspace_save(solve,pp,denom)

    type(soltyp), intent(inout) :: solve
    real(rp),     intent(in)    :: pp(*)
    real(rp),     intent(in)    :: denom
    integer(ip)                 :: ii,kk

    if( solve % kfl_save_krylov > 0 ) then

       solve % krylov_size = solve % krylov_size + 1
       if( solve % krylov_size > solve % kfl_save_krylov ) then
          if( INOTMASTER ) then
             do kk = 1,solve % kfl_save_krylov-1
                do ii = 1,solve % nunkn
                   solve % krylov_dir(ii,kk) = solve % krylov_dir(ii,kk+1)
                end do
             end do
          end if
          solve % krylov_size = solve % kfl_save_krylov
       end if

       if( INOTMASTER ) then
          kk = solve % krylov_size
          do ii = 1,solve % nunkn
             solve % krylov_dir(ii,kk) = pp(ii)    ! pk
          end do
          solve % krylov_dot_products(kk) = denom  ! (Apk,pk)
       end if

    end if

  end subroutine solver_krylov_subspace_save

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    18/05/2017
  !> @brief   Initial guess using a previous a Krylov subspace
  !> @details Project initial solution on last Krylov supbspace to compute
  !>          an initial guess to solve Ax=b
  !>          Demonstration:
  !>          \verbatim
  !>
  !>          To take into account initial solution:
  !>          Let
  !>          (1) z = x - x0
  !>          so that z0 = 0 Ax = b => A (z+x0) = b =>
  !>          (2) A z = b0 is the new system
  !>          (3) b0 = b - Ax0 = r0
  !>
  !>          We have
  !>
  !>          z        = sum_i ai pi
  !>
  !>          A.z      = sum_i ai A.pi       (multiply by A)
  !>          (pk,A.z) = sum_i ai (pk,A.pi)  (dot product with pk)
  !>
  !>          Vectors p's are conjugate so (pk,A.pi) = (pk,A.pk) \delta_ik
  !>
  !>          (pk,A.z) = ak (pk,A.pk)
  !>          (pk,b0)  = ak (pk,A.pk)        (Equation 2)
  !>          (pk,r0)  = ak (pk,A.pk)        (Equation 3)
  !>
  !>          and therefore:
  !>
  !>                 (pk,r0)
  !>          ak =  ---------
  !>                (pk,A.pk)
  !>
  !>          Using Equation 1 , we thus have:
  !>
  !>                           (pk,r0)
  !>          x = x0 + sum_k  --------- pk
  !>                          (pk,A.pk)
  !>
  !>          \endverbatim
  !>
  !>          References:
  !>          [1] https://en.wikipedia.org/wiki/Conjugate_gradient_method
  !>          [2] D. Tromeur-Dervout, Y. Vassilevski, Choice of initial
  !>              guess in iterative solution of series of systems arising
  !>              in fluid flow simulations,
  !>              J. Compt. Phys. 219, 210-227 (2006)
  !>
  !>
  !----------------------------------------------------------------------

  subroutine solver_krylov_subspace_initial_guess(solve,xx,rr)

    type(soltyp), intent(inout) :: solve
    real(rp),     intent(inout) :: xx(*)
    real(rp),     intent(in)    :: rr(*)
    integer(ip)                 :: ii,kk,nbnodes,nbvar,nrows
    real(rp)                    :: alpha,rho,dumm1(2),dumm2(2)

    if( solve % kfl_save_krylov > 0 ) then

       nbnodes = solve % nequa
       nbvar   = solve % ndofn
       nrows   = nbnodes * nbvar

       if( INOTMASTER ) then

          do kk = 1,solve % krylov_size
             call solver_parallel_scalar_product(solve,rr,solve % krylov_dir(:,kk),rho) ! alpha = (z,p)
             alpha = rho / solve % krylov_dot_products(kk)                              ! alpha = (z,p)/(Ap,p)
             do ii = 1,nrows
                xx(ii) = xx(ii) + alpha * solve % krylov_dir(ii,kk)                     ! x = x + alphai pi
             end do
          end do

       else

          do kk = 1,solve % krylov_size             
             call solver_parallel_scalar_product(solve,dumm1,dumm2,rho)
          end do

       end if

       solve % krylov_size = 0

    end if

  end subroutine solver_krylov_subspace_initial_guess

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    21/12/2016
  !> @brief   Solve a system with mass matrix
  !> @details Solve a diagonal system using the lumped mass matrix:
  !>
  !>          M U = B
  !>
  !>          RHS is U
  !>          Solution is stored in U
  !>
  !----------------------------------------------------------------------

  subroutine solver_lumped_mass_system(ndofn,u,EXCHANGE)

    implicit none
    integer(ip), intent(in)              :: ndofn              !< run through either npoin or nbopo
    real(rp),    intent(inout), target   :: u(ndofn,*)         !< in-out corrected vector
    logical(lg), intent(in),    optional :: EXCHANGE
    integer(ip)                          :: ii,icoup
    real(rp)                             :: dummr,u1(2),u2(2)
    logical(lg)                          :: if_exchange
    
    if( INOTMASTER ) then
       !
       ! Check if we should exchange
       !
       if_exchange = .true.
       if( present(EXCHANGE) ) if_exchange = EXCHANGE

       if( if_exchange ) then
          !
          ! Periodicity
          !
          call solver_periodicity_rhs(ndofn,u)
          !
          ! Parall service: exchange result between slaves
          !
          call pararr('SLX',NPOIN_TYPE,npoin*ndofn,u)
          !call PAR_INTERFACE_NODE_EXCHANGE(ndofn,u,'SUM','IN MY CODE')
          !
          ! Periodicity: recover solution on slave nodes
          !
          call solver_periodicity_solution(ndofn,u)
       end if
       !
       ! Solve system
       !
       if(      ndofn == 1 ) then
          do ii = 1,npoin
             u(1,ii) = u(1,ii) / vmass(ii)
          end do
       else if( ndofn == 2 ) then
          do ii = 1,npoin
             u(1,ii) = u(1,ii) / vmass(ii)
             u(2,ii) = u(2,ii) / vmass(ii)
          end do
       else if( ndofn == 3 ) then
          do ii = 1,npoin
             u(1,ii) = u(1,ii) / vmass(ii)
             u(2,ii) = u(2,ii) / vmass(ii)
             u(3,ii) = u(3,ii) / vmass(ii)
          end do
       else
          do ii = 1,npoin
             u(1:ndofn,ii) = u(1:ndofn,ii) / vmass(ii)
          end do
       end if
    end if
    ! 
    ! Impose Dirichlet condition
    !
    if( ncoup_implicit_d > 0 ) then
       do ii = 1,ncoup_implicit_d
          icoup = lcoup_implicit_d(ii) 
          if( I_AM_IN_COUPLING(icoup) ) then
             if( IMASTER ) then
                call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,u1)
             else
                call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,u)
             end if
          end if
       end do
    end if

  end subroutine solver_lumped_mass_system

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    21/12/2016
  !> @brief   Solve a diagonal system
  !> @details Solve a diagonal system:
  !>
  !>          D U = B
  !>
  !>          RHS is U
  !>          Solution is stored in U
  !>
  !----------------------------------------------------------------------

  subroutine solver_diagonal_system(ndofn,diag,u)

    implicit none
    integer(ip), intent(in)            :: ndofn            !< run through either npoin or nbopo
    real(rp),    intent(in)            :: diag(*)          !< in-out corrected vector
    real(rp),    intent(inout), target :: u(ndofn,*)       !< in-out corrected vector
    integer(ip)                        :: ii,icoup
    real(rp)                           :: dummr,u1(2),u2(2)
    real(rp),    pointer               :: u_cpy(:,:)
    
    if( INOTMASTER ) then
       !
       ! Periodicity
       !
       call solver_periodicity_rhs(ndofn,u)
       !
       ! Parall service: exchange result between slaves
       !
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,u,'SUM','IN MY CODE')
       !
       ! Periodicity: recover solution on slave nodes
       !
       call solver_periodicity_solution(ndofn,u)
       !
       ! Solve system
       !
       do ii = 1,npoin
          u(1:ndofn,ii) = u(1:ndofn,ii)  / diag(ii)
       end do
    end if
    
    if( ncoup_implicit_d > 0 ) then
       if( IMASTER ) then
          do ii = 1,ncoup_implicit_d
             icoup = lcoup_implicit_d(ii) 
             if( I_AM_IN_COUPLING(icoup) ) then
                call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,u1,u2)
             end if
          end do
       else
          nullify(u_cpy)
          allocate(u_cpy(ndofn,npoin)) 
          u_cpy(1:ndofn,1:npoin) = u(1:ndofn,1:npoin)
          do ii = 1,ncoup_implicit_d
             icoup = lcoup_implicit_d(ii) 
             if( I_AM_IN_COUPLING(icoup) ) then
                call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,u,u_cpy)
             end if
          end do
          deallocate(u_cpy)
       end if
    end if


  end subroutine solver_diagonal_system

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    21/12/2016
  !> @brief   Solve a system with mass matrix
  !> @details Solve a diagonal system using the close mass matrix:
  !>
  !>          M U = B
  !>
  !>          RHS is U
  !>          Solution is stored in U
  !>
  !----------------------------------------------------------------------

  subroutine solver_close_mass_system(ndofn,u)

    implicit none
    integer(ip), intent(in)            :: ndofn              !< run through either npoin or nbopo
    real(rp),    intent(inout), target :: u(ndofn,*)         !< in-out corrected vector
    integer(ip)                        :: ii
    real(rp)                           :: dummr

    if( INOTMASTER ) then
       !
       ! Periodicity
       !
       call solver_periodicity_rhs(ndofn,u)
       !
       ! Parall service: exchange result between slaves
       !
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,u,'SUM','IN MY CODE')
       !
       ! Periodicity: recover solution on slave nodes
       !
       call solver_periodicity_solution(ndofn,u)
       !
       ! Solve system
       !
       do ii = 1,npoin
          u(1:ndofn,ii) = u(1:ndofn,ii)  / vmasc(ii)
       end do

    end if

  end subroutine solver_close_mass_system

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    11/12/2017
  !> @brief   Residual and residual norm
  !> @details Compute the residual of an equation:
  !>          r = b - Ax
  !>          It is assumed that b was already assembled between
  !>          neighbors
  !>
  !----------------------------------------------------------------------

  subroutine solver_parallel_residual_and_norm(&
       solve,an,xx,bb,RESIDUAL,RESIDUAL_NORM)
    type(soltyp), intent(inout)                    :: solve          !< Solver structure
    real(rp),     intent(in)                       :: an(*)          !< Matrix
    real(rp),     intent(inout)                    :: xx(*)          !< Solution
    real(rp),     intent(in)                       :: bb(*)          !< Solution
    real(rp),     intent(inout), pointer, optional :: RESIDUAL(:)    !< Residual
    real(rp),     intent(out),            optional :: RESIDUAL_NORM  !< Residual norm
    integer(ip)                                    :: ndof,nn,nsize
    real(rp),     pointer                          :: rr(:)
     
    nn    = solve % nequa
    ndof  = solve % ndofn
    nsize = nn * ndof

    if( present(RESIDUAL) ) then
       if( INOTMASTER ) then
          if( size(RESIDUAL,KIND=ip) < nsize ) call runend('SOLVER_PARALLEL_RESIDUAL_AND_NORM: WRONG RESIDUAL')
          call solver_parallel_SpMV(solve,an,xx,RESIDUAL,OPENMP=.true.)
          RESIDUAL(1:nsize) = bb(1:nsize)-RESIDUAL(1:nsize)
       end if
       if( present(RESIDUAL_NORM) ) call solver_parallel_vector_L2norm(solve,RESIDUAL,RESIDUAL_NORM,OPENMP=.true.) 
    else
       allocate(rr(max(1_ip,nsize)))  
       if( INOTMASTER ) then
          call solver_parallel_SpMV(solve,an,xx,rr,OPENMP=.true.)
          rr(1:nsize) = bb(1:nsize)-rr(1:nsize)
       end if
       if( present(RESIDUAL_NORM) ) call solver_parallel_vector_L2norm(solve,rr,RESIDUAL_NORM,OPENMP=.true.)
       deallocate(rr)      
       
     end if
    
   end subroutine solver_parallel_residual_and_norm
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    00/05/2016
  !> @brief   Matrix vector product
  !> @details Compute the residual norm of an equation
  !>
  !----------------------------------------------------------------------

  subroutine solver_array_residual_norm(&
       solve,norm,nn1,nn2,v1,v2,kcom1,kcom2,kdime,relax,redif)

    type(soltyp), intent(in)  :: solve
    integer(ip),  intent(in)  :: nn1
    integer(ip),  intent(in)  :: nn2
    integer(ip),  intent(in)  :: norm
    integer(ip),  intent(in)  :: kcom1
    integer(ip),  intent(in)  :: kcom2
    integer(ip),  intent(in)  :: kdime
    real(rp),     intent(in)  :: v1(*)
    real(rp),     intent(in)  :: v2(*)
    real(rp),     intent(in)  :: relax
    real(rp),     intent(out) :: redif
    integer(ip)               :: ii,i1,i2,idime
    real(rp)                  :: va,vo,vd,ra,ro
    real(rp)                  :: resid(2)

    ra    = relax
    ro    = 1.0_rp-ra
    resid = 0.0_rp

    select case ( norm )

    case ( 0_ip )
       !
       ! Linf norm
       !
       if( INOTMASTER ) then
          do ii = 1,solve % nequ1
             i1 = (ii-1) * nn1 + kcom1
             i2 = (ii-1) * nn2 + kcom2
             do idime = 0,kdime-1
                vo       = v2(i2+idime)
                va       = ra*v1(i1+idime) + ro*vo
                resid(1) = max(resid(1),abs(va-vo))
                resid(2) = max(resid(2),abs(va))
             end do
          end do
          do ii = solve % nequ2,solve % nequ3
             i1 = (ii-1) * nn1 + kcom1
             i2 = (ii-1) * nn2 + kcom2
             do idime = 0,kdime-1
                vo       = v2(i2+idime)
                va       = ra*v1(i1+idime)+ro*vo
                resid(1) = max(resid(1),abs(va-vo))
                resid(2) = max(resid(2),abs(va))
             end do
          end do
       end if

       call PAR_MAX(2_ip,resid)

       if( resid(2) > zeror ) then
          redif = resid(1) / resid(2)
       else
          redif = resid(1)
       end if

    case ( 1_ip )
       !
       ! L1 norm
       !
       if( INOTMASTER ) then
          do ii = 1,solve % nequ1
             i1 = (ii-1) * nn1 + kcom1
             i2 = (ii-1) * nn2 + kcom2
             do idime = 0,kdime-1
                vo       = v2(i2+idime)
                va       = ra*v1(i1+idime)+ro*vo
                resid(1) = resid(1) + abs(va-vo)
                resid(2) = resid(2) + abs(va)
             end do
          end do
          do ii = solve % nequ2,solve % nequ3
             i1 = (ii-1) * nn1 + kcom1
             i2 = (ii-1) * nn2 + kcom2
             do idime = 0,kdime-1
                vo       = v2(i2+idime)
                va       = ra*v1(i1+idime)+ro*vo
                resid(1) = resid(1) + abs(va-vo)
                resid(2) = resid(2) + abs(va)
             end do
          end do
       end if

       call PAR_SUM(2_ip,resid)

       if( resid(2) > zeror ) then
          redif = resid(1) / resid(2)
       else
          redif = resid(1)
       end if

    case ( 2_ip )
       !
       ! L2 norm
       !
       if( INOTMASTER ) then
          do ii = 1,solve % nequ1
             i1 = (ii-1) * nn1 + kcom1
             i2 = (ii-1) * nn2 + kcom2
             do idime = 0,kdime-1
                vo       = v2(i2+idime)
                va       = ra*v1(i1+idime)+ro*vo
                vd       = va-vo
                resid(1) = resid(1) + vd*vd
                resid(2) = resid(2) + va*va
             end do
          end do
          do ii = solve % nequ2,solve % nequ3
             i1 = (ii-1) * nn1 + kcom1
             i2 = (ii-1) * nn2 + kcom2
             do idime = 0,kdime-1
                vo       = v2(i2+idime)
                va       = ra*v1(i1+idime)+ro*vo
                vd       = va-vo
                resid(1) = resid(1) + vd*vd
                resid(2) = resid(2) + va*va
             end do
          end do
       end if

       call PAR_SUM(2_ip,resid)

       if( resid(2) > zeror ) then
          redif = sqrt(resid(1)/resid(2))
       else
          redif = sqrt(resid(1))
       end if

    end select

  end subroutine solver_array_residual_norm

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2017
  !> @brief   L2 norm of a vector
  !> @details Compute the L2 norm XNORM of a vector XX
  !>
  !----------------------------------------------------------------------

  subroutine solver_parallel_vector_L2norm_1(solve,xx,xnorm,MY_TIMING,OPENMP)
    type(soltyp), intent(inout)         :: solve
    real(rp),     intent(in)            :: xx(*)
    real(rp),     intent(out)           :: xnorm
    real(rp),     intent(out), optional :: MY_TIMING(2)
    logical(lg),  intent(in),  optional :: OPENMP
    
    call solver_parallel_scalar_product(solve,xx,xx,xnorm,MY_TIMING,OPENMP)
    xnorm = sqrt(xnorm)
    
  end subroutine solver_parallel_vector_L2norm_1
  
  subroutine solver_parallel_vector_L2norm_2(solve,xx,yy,xnorm,ynorm,MY_TIMING,OPENMP)
    type(soltyp), intent(inout)         :: solve
    real(rp),     intent(in)            :: xx(*)
    real(rp),     intent(in)            :: yy(*)
    real(rp),     intent(out)           :: xnorm
    real(rp),     intent(out)           :: ynorm
    real(rp),     intent(out), optional :: MY_TIMING(2)
    logical(lg),  intent(in),  optional :: OPENMP
    
    call solver_parallel_scalar_product(solve,xx,xx,yy,yy,xnorm,ynorm,MY_TIMING,OPENMP)
    xnorm = sqrt(xnorm)
    ynorm = sqrt(ynorm)
    
  end subroutine solver_parallel_vector_L2norm_2
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2017
  !> @brief   Linf norm of a vector
  !> @details Compute the Linf norm XNORM of a vector XX
  !>
  !----------------------------------------------------------------------

  subroutine solver_parallel_vector_Linfnorm(solve,xx,xnorm)
    type(soltyp), intent(in)            :: solve
    real(rp),     intent(in),  pointer  :: xx(:)
    real(rp),     intent(out)           :: xnorm
    integer(ip)                         :: nsize
    
    nsize = solve % nequa_own * solve % ndofn
    xnorm = maxval(xx(1:nsize))
    call PAR_MAX(xnorm,'IN MY CODE')

  end subroutine solver_parallel_vector_Linfnorm
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Scalar product
  !> @details Compute the parallel scalar product XXDOTYY = XX.YY
  !>
  !----------------------------------------------------------------------

  subroutine solver_parallel_scalar_product_t(solve,xx,yy,xxdotyy,MY_TIMING,OPENMP,MPI)
    type(soltyp), intent(inout)         :: solve
    real(rp),     intent(in)            :: xx(*)
    real(rp),     intent(in)            :: yy(*)
    real(rp),     intent(out)           :: xxdotyy
    real(rp),     intent(out), optional :: MY_TIMING(2)
    logical(lg),  intent(in),  optional :: OPENMP
    logical(lg),  intent(in),  optional :: MPI
    integer(ip)                         :: ii,jj,kk,nsize
    real(rp)                            :: time1,time2,time3
    logical(lg)                         :: use_openmp
    logical(lg)                         :: use_mpi

#ifdef BLAS
    real(rp) :: DDOT
    external DDOT
#endif
    call cputim(time1)
    
    nsize      = solve % nequa_own * solve % ndofn
    use_openmp = .false.
    use_mpi    = .true.
    xxdotyy    = 0.0_rp
    
    if( present(OPENMP) ) use_openmp = OPENMP
    if( present(MPI) )    use_mpi    = MPI

    if( nsize > 0 .and. INOTMASTER ) then

       if( solve % kfl_mask == 0 ) then
#ifdef BLAS
          xxdotyy = DDOT(nsize,xx,1_ip,yy,1_ip)
#else
          if( use_openmp ) then
             xxdotyy = 0.0_rp
             !$OMP PARALLEL    DO SCHEDULE (STATIC) &
             !$OMP DEFAULT   ( NONE )               &
             !$OMP PRIVATE   ( ii )                 &
             !$OMP SHARED    ( xx, yy, nsize )      &
             !$OMP REDUCTION ( +:xxdotyy )  
             do ii = 1,nsize
                xxdotyy = xxdotyy + xx(ii) * yy(ii)
             end do
             !$OMP END PARALLEL DO
          else
             xxdotyy = dot_product(xx(1:nsize),yy(1:nsize))
          end if
#endif
       else

          xxdotyy = 0.0_rp
          kk = 0
          do ii = 1,solve % nequa_own
             do jj = 1,solve % ndofn
                kk = kk + 1
                xxdotyy = xxdotyy + xx(kk) * yy(kk) * solve % mask(ii)
             end do
          end do
          
       end if
       
    end if
    
    call cputim(time2)
    if( use_mpi ) call PAR_SUM(xxdotyy,'IN MY CODE')

    call cputim(time3)
    if( present(MY_TIMING) ) then
       MY_TIMING(1) = time2 - time1
       MY_TIMING(2) = time3 - time2
    else
       solve % cpu_dot(1) = solve % cpu_dot(1) + time2 - time1
       solve % cpu_dot(2) = solve % cpu_dot(2) + time3 - time2
       solve % cpu_dot(3) = solve % cpu_dot(1) + solve % cpu_dot(2)
       solve % num_dot    = solve % num_dot + 1
    end if

  end subroutine solver_parallel_scalar_product_t

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Double scalar product
  !> @details Compute XX.YY and XX.ZZ
  !>
  !----------------------------------------------------------------------

  subroutine solver_parallel_double_scalar_product(solve,xx,yy,zz,xxdotyy,xxdotzz,MY_TIMING,OPENMP)
    type(soltyp), intent(inout)         :: solve
    real(rp),     intent(in)            :: xx(*)
    real(rp),     intent(in)            :: yy(*)
    real(rp),     intent(in)            :: zz(*)
    real(rp),     intent(out)           :: xxdotyy
    real(rp),     intent(out)           :: xxdotzz
    real(rp),     intent(out), optional :: MY_TIMING(2)
    logical(lg),  intent(in),  optional :: OPENMP
    integer(ip)                         :: ii,nsize
    real(rp)                            :: time1,time2,time3
    logical(lg)                         :: use_openmp
    real(rp)                            :: dots(2)

#ifdef BLAS
    real(rp) :: DDOT
    external DDOT
#endif
    call cputim(time1)
    
    nsize      = solve % nequa_own * solve % ndofn
    use_openmp = .false.
    xxdotyy    = 0.0_rp
    xxdotzz    = 0.0_rp
    
    if( present(OPENMP) ) use_openmp = OPENMP

    if( nsize > 0 .and. INOTMASTER ) then

#ifdef BLAS
       xxdotyy = DDOT(nsize,xx,1_ip,yy,1_ip)
       xxdotzz = DDOT(nsize,xx,1_ip,zz,1_ip)
#else
       if( use_openmp ) then
          xxdotyy = 0.0_rp
          xxdotzz = 0.0_rp
          !$OMP PARALLEL    DO SCHEDULE (STATIC) &
          !$OMP DEFAULT   ( NONE )               &
          !$OMP PRIVATE   ( ii )                 &
          !$OMP SHARED    ( xx, yy, zz, nsize )  &
          !$OMP REDUCTION ( +:xxdotyy,xxdotzz )  
          do ii = 1,nsize
             xxdotyy = xxdotyy + xx(ii) * yy(ii)
             xxdotzz = xxdotzz + xx(ii) * zz(ii)
          end do
          !$OMP END PARALLEL DO
       else
          xxdotyy = dot_product(xx(1:nsize),yy(1:nsize))
          xxdotzz = dot_product(xx(1:nsize),zz(1:nsize))
       end if
#endif
    end if

    call cputim(time2)
    dots(1) = xxdotyy
    dots(2) = xxdotzz
    call PAR_SUM(2_ip,dots)
    xxdotyy = dots(1) 
    xxdotzz = dots(2) 

    call cputim(time3)
    if( present(MY_TIMING) ) then
       MY_TIMING(1) = time2 - time1
       MY_TIMING(2) = time3 - time2
    else
       solve % cpu_dot(1) = solve % cpu_dot(1) + time2 - time1
       solve % cpu_dot(2) = solve % cpu_dot(2) + time3 - time2
       solve % cpu_dot(3) = solve % cpu_dot(1) + solve % cpu_dot(2)
       solve % num_dot    = solve % num_dot + 1
    end if

  end subroutine solver_parallel_double_scalar_product

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Double scalar product
  !> @details Compute WW.XX and YY.ZZ
  !>
  !----------------------------------------------------------------------

  subroutine solver_parallel_double_scalar_product_2(solve,ww,xx,yy,zz,wwdotxx,yydotzz,MY_TIMING,OPENMP)
    type(soltyp), intent(inout)         :: solve
    real(rp),     intent(in)            :: ww(*)
    real(rp),     intent(in)            :: xx(*)
    real(rp),     intent(in)            :: yy(*)
    real(rp),     intent(in)            :: zz(*)
    real(rp),     intent(out)           :: wwdotxx
    real(rp),     intent(out)           :: yydotzz
    real(rp),     intent(out), optional :: MY_TIMING(2)
    logical(lg),  intent(in),  optional :: OPENMP
    integer(ip)                         :: ii,nsize
    real(rp)                            :: time1,time2,time3
    logical(lg)                         :: use_openmp
    real(rp)                            :: dots(2)

#ifdef BLAS
    real(rp) :: DDOT
    external DDOT
#endif
    call cputim(time1)
    nsize = solve % nequa_own * solve % ndofn
    use_openmp = .false.
    if( present(OPENMP) ) use_openmp = OPENMP

    if( nsize > 0 .and. INOTMASTER ) then

#ifdef BLAS
       wwdotxx = DDOT(nsize,ww,1_ip,xx,1_ip)
       yydotzz = DDOT(nsize,yy,1_ip,zz,1_ip)
#else
       if( use_openmp ) then
          wwdotxx = 0.0_rp
          yydotzz = 0.0_rp
          !$OMP PARALLEL    DO SCHEDULE (STATIC)     &
          !$OMP DEFAULT   ( NONE )                   &
          !$OMP PRIVATE   ( ii )                     &
          !$OMP SHARED    ( ww, xx, yy, zz, nsize )  &
          !$OMP REDUCTION ( +:wwdotxx,yydotzz )  
          do ii = 1,nsize
             wwdotxx = wwdotxx + ww(ii) * xx(ii)
             yydotzz = yydotzz + yy(ii) * zz(ii)
          end do
          !$OMP END PARALLEL DO
       else
          wwdotxx = dot_product(ww(1:nsize),xx(1:nsize))
          yydotzz = dot_product(yy(1:nsize),zz(1:nsize))
       end if
#endif
    end if

    call cputim(time2)
    dots(1) = wwdotxx
    dots(2) = yydotzz
    call PAR_SUM(2_ip,dots)
    wwdotxx = dots(1) 
    yydotzz = dots(2) 

    call cputim(time3)
    if( present(MY_TIMING) ) then
       MY_TIMING(1) = time2 - time1
       MY_TIMING(2) = time3 - time2
    else
       solve % cpu_dot(1) = solve % cpu_dot(1) + time2 - time1
       solve % cpu_dot(2) = solve % cpu_dot(2) + time3 - time2
       solve % cpu_dot(3) = solve % cpu_dot(1) + solve % cpu_dot(2)
       solve % num_dot    = solve % num_dot + 1
    end if

  end subroutine solver_parallel_double_scalar_product_2

  subroutine solver_parallel_scalar_product_1(nbvar,xx,yy,xxdotyy)
    integer(ip), intent(in)           :: nbvar
    real(rp),    intent(in),  pointer :: xx(:)
    real(rp),    intent(in),  pointer :: yy(:)
    real(rp),    intent(out)          :: xxdotyy
    integer(ip)                       :: ii,nsize
#ifdef BLAS
    real(rp) :: DDOT
    external DDOT
#endif

    nsize = nbvar * npoi3
    
    if( INOTMASTER .and. nsize > 0 ) THEN
#ifdef BLAS
       xxdotyy = DDOT(nsize,xx,1_ip,yy,1_ip)
#else
       xxdotyy = 0.0_rp
       do ii = 1,nsize
          xxdotyy = xxdotyy + xx(ii) * yy(ii)
       end do
#endif       
    end if

    call PAR_SUM(xxdotyy,'IN MY CODE')

  end subroutine solver_parallel_scalar_product_1

  subroutine solver_parallel_scalar_product_2_3(nbvar,xx,yy,zz,xxdotyyzz,wsynch)
    integer(ip), intent(in)            :: nbvar
    real(rp),    intent(in),  pointer  :: xx(:)
    real(rp),    intent(in),  pointer  :: yy(:)
    real(rp),    intent(in),  pointer  :: zz(:)
    real(rp),    intent(out)           :: xxdotyyzz(2)
    character(*),intent(in),  optional :: wsynch
    integer(ip)                        :: ii
    logical(lg)                        :: asynch
    logical(lg)                        :: all_reduce
#ifdef BLAS
    real(rp) :: DDOT
    external DDOT
#endif

    xxdotyyzz  = 0.0_rp
    asynch     = .false.
    all_reduce = .true.
    if( present(wsynch) ) then
       if( trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
       else if( trim(wsynch) == 'DO NOT REDUCE' ) then
          all_reduce = .false.
       end if
    end if

    if( INOTMASTER ) then
       !
       ! Loop over interior nodes
       !
#ifdef BLAS
       xxdotyyzz(1) = DDOT(nbvar*npoi1,xx,1_ip,yy,1_ip)
       xxdotyyzz(2) = DDOT(nbvar*npoi1,xx,1_ip,zz,1_ip)
#else
       do ii = 1,nbvar*npoi1
          xxdotyyzz(1) = xxdotyyzz(1) + xx(ii) * yy(ii)
          xxdotyyzz(2) = xxdotyyzz(2) + xx(ii) * zz(ii)
       end do
#endif
       !
       ! Loop over own boundary nodes
       !
       do ii = (npoi2-1)*nbvar+1,npoi3*nbvar
          xxdotyyzz(1) = xxdotyyzz(1) + xx(ii) * yy(ii)
          xxdotyyzz(2) = xxdotyyzz(2) + xx(ii) * zz(ii)
       end do
    end if

    if( all_reduce .and. IPARALL ) then
       if( asynch ) then
          call PAR_SUM(2_ip,xxdotyyzz,'IN MY CODE',wsynch)
       else
          call PAR_SUM(2_ip,xxdotyyzz,'IN MY CODE')
       end if
    end if

  end subroutine solver_parallel_scalar_product_2_3

  subroutine solver_parallel_scalar_product_2_4(nbvar,ww,xx,yy,zz,wwxx_yyzz,wsynch)
    integer(ip), intent(in)            :: nbvar
    real(rp),    intent(in),  pointer  :: ww(:)
    real(rp),    intent(in),  pointer  :: xx(:)
    real(rp),    intent(in),  pointer  :: yy(:)
    real(rp),    intent(in),  pointer  :: zz(:)
    real(rp),    intent(out)           :: wwxx_yyzz(2)
    character(*),intent(in),  optional :: wsynch
    integer(ip)                        :: ii
    logical(lg)                        :: asynch
    logical(lg)                        :: all_reduce

    wwxx_yyzz  = 0.0_rp
    asynch     = .false.
    all_reduce = .true.
    if( present(wsynch) ) then
       if( trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
       else if( trim(wsynch) == 'DO NOT REDUCE' ) then
          all_reduce = .false.
       end if
    end if

    if( INOTMASTER ) then
       !
       ! Loop over interior nodes
       !
       do ii = 1,nbvar*npoi1
          wwxx_yyzz(1) = wwxx_yyzz(1) + ww(ii) * xx(ii)
          wwxx_yyzz(2) = wwxx_yyzz(2) + yy(ii) * zz(ii)
       end do
       !
       ! Loop over own boundary nodes
       !
       do ii = (npoi2-1)*nbvar+1,npoi3*nbvar
          wwxx_yyzz(1) = wwxx_yyzz(1) + ww(ii) * xx(ii)
          wwxx_yyzz(2) = wwxx_yyzz(2) + yy(ii) * zz(ii)
       end do
    end if

    if( all_reduce .and. IPARALL ) then
       if( asynch ) then
          call PAR_SUM(2_ip,wwxx_yyzz,'IN MY CODE',wsynch)
       else
          call PAR_SUM(2_ip,wwxx_yyzz,'IN MY CODE')
       end if
    end if

  end subroutine solver_parallel_scalar_product_2_4

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Matrix vector product
  !> @details Multiply a non symmetric matrix stored by a vector
  !>
  !----------------------------------------------------------------------

  subroutine solver_smvp(solve,an,xx,yy)
    type(soltyp), intent(inout)  :: solve
    real(rp),     intent(in)     :: an(*)
    real(rp),     intent(in)     :: xx(*)
    real(rp),     intent(out)    :: yy(*)

    if(      solve % kfl_format == SOL_BCSR_FORMAT ) then
       !
       ! BCSR format
       !
       call solver_CSR_smvp(solve % nequa,solve % ndofn,solve % ndofn,an,solve % ja,solve % ia,xx,yy)

    else if( solve % kfl_format == SOL_COO_FORMAT ) then
       !
       ! COO format
       !
       call solver_COO_smvp(solve % nequa,solve % nnz,solve % ndofn,solve % ndofn,an,solve % rows,solve % cols,xx,yy)

    end if

  end subroutine solver_smvp

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Matrix vector product
  !> @details Multiply a non symmetric matrix stored in BCSR by a vector
  !>          YY = A XX
  !>          INPUT
  !>          NBNODES .... Number of equations
  !>            NBVAR ...... Number of variables
  !>            AN ......... Matrix
  !>            JA ......... List of elements
  !>            IA ......... Pointer to list of elements
  !>            XX ......... Vector
  !>          OUTPUT
  !>            YY ......... result vector
  !>
  !----------------------------------------------------------------------

  subroutine solver_CSR_smvp(nbnodes,ndof1,ndof2,an,ja,ia,xx,yy)
    integer(ip),  intent(in)           :: nbnodes
    integer(ip),  intent(in)           :: ndof1
    integer(ip),  intent(in)           :: ndof2
    real(rp),     intent(in)           :: an(ndof1,ndof2,*)
    integer(ip),  intent(in)           :: ja(*),ia(*)
    real(rp),     intent(in)           :: xx(ndof1,*)
    real(rp),     intent(out)          :: yy(ndof2,*)
    integer(ip)                        :: ii,jj,kk,ll,col
    real(rp)                           :: raux,raux1,raux2,raux3

    if( INOTMASTER ) then

       if( ndof1 == 1 .and. ndof2 == 1 ) then
          !
          ! NBVAR=1
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)  &
          !$OMP DEFAULT  ( NONE )                &
          !$OMP PRIVATE  ( ii, jj, col, raux )   &
          !$OMP SHARED   ( xx, yy, an, nbnodes )
          do ii = 1,nbnodes
             yy(1,ii) = 0.0_rp
             do jj   = ia(ii),ia(ii+1)-1
                col  = ja(jj)
                raux = xx(1,col)
                yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
             end do
          end do
          !$OMP END PARALLEL DO

       else if( ndof1 == 2 .and. ndof2 == 2 ) then
          !
          ! NBVAR=2
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)        &
          !$OMP DEFAULT  ( NONE )                      &
          !$OMP PRIVATE  ( ii, jj, col, raux1, raux2 ) &
          !$OMP SHARED   ( xx, yy, an, nbnodes )
          do ii = 1,nbnodes
             yy(1,ii) = 0.0_rp
             yy(2,ii) = 0.0_rp
             do jj       = ia(ii),ia(ii+1)-1
                col      = ja(jj)
                raux1    = xx(1,col)
                raux2    = xx(2,col)
                yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
                yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
                yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
                yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
             end do
          end do
          !$OMP END PARALLEL DO

       else if( ndof1 == 3 .and. ndof2 == 3 ) then
          !
          ! NBVAR=3
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)               &
          !$OMP DEFAULT  ( NONE )                             &
          !$OMP PRIVATE  ( ii, jj, col, raux1, raux2, raux3 ) &
          !$OMP SHARED   ( xx, yy, an, nbnodes )
          do ii = 1,nbnodes
             yy(1,ii) = 0.0_rp
             yy(2,ii) = 0.0_rp
             yy(3,ii) = 0.0_rp
             do jj       = ia(ii),ia(ii+1)-1
                col      = ja(jj)
                raux1    = xx(1,col)
                raux2    = xx(2,col)
                raux3    = xx(3,col)
                yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux1
                yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux2
                yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux3
                yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux1
                yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux2
                yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux3
                yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux1
                yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux2
                yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux3
             end do
          end do
          !$OMP END PARALLEL DO

       else if( ndof1 == 4 .and. ndof2 == 4 ) then
          !
          ! NBVAR=4
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)  &
          !$OMP DEFAULT  ( NONE )                &
          !$OMP PRIVATE  ( ii, jj, col, raux )   &
          !$OMP SHARED   ( xx, yy, an, nbnodes )
          do ii = 1,nbnodes
             yy(1,ii) = 0.0_rp
             yy(2,ii) = 0.0_rp
             yy(3,ii) = 0.0_rp
             yy(4,ii) = 0.0_rp
             do jj       = ia(ii),ia(ii+1)-1
                col      = ja(jj)
                raux     = xx(1,col)
                yy(1,ii) = yy(1,ii) + an(1,1,jj) * raux
                yy(2,ii) = yy(2,ii) + an(1,2,jj) * raux
                yy(3,ii) = yy(3,ii) + an(1,3,jj) * raux
                yy(4,ii) = yy(4,ii) + an(1,4,jj) * raux
                raux     = xx(2,col)
                yy(1,ii) = yy(1,ii) + an(2,1,jj) * raux
                yy(2,ii) = yy(2,ii) + an(2,2,jj) * raux
                yy(3,ii) = yy(3,ii) + an(2,3,jj) * raux
                yy(4,ii) = yy(4,ii) + an(2,4,jj) * raux
                raux     = xx(3,col)
                yy(1,ii) = yy(1,ii) + an(3,1,jj) * raux
                yy(2,ii) = yy(2,ii) + an(3,2,jj) * raux
                yy(3,ii) = yy(3,ii) + an(3,3,jj) * raux
                yy(4,ii) = yy(4,ii) + an(3,4,jj) * raux
                raux     = xx(4,col)
                yy(1,ii) = yy(1,ii) + an(4,1,jj) * raux
                yy(2,ii) = yy(2,ii) + an(4,2,jj) * raux
                yy(3,ii) = yy(3,ii) + an(4,3,jj) * raux
                yy(4,ii) = yy(4,ii) + an(4,4,jj) * raux
             end do
          end do
          !$OMP END PARALLEL DO

       else
          !
          ! NBVAR = whatever
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)              &
          !$OMP DEFAULT  ( NONE )                            &
          !$OMP PRIVATE  ( ii, jj, ll, col, raux )           &
          !$OMP SHARED   ( xx, yy, an, ndof1, ndof2, nbnodes )
          do ii = 1,nbnodes
             yy(1:ndof2,ii) = 0.0_rp
             do jj  = ia(ii),ia(ii+1)-1
                col = ja(jj)
                do ll = 1,ndof1
                   raux = xx(ll,col)
                   yy(1:ndof2,ii) = yy(1:ndof2,ii) + an(ll,1:ndof2,jj) * raux
                end do
             end do
          end do
          !$OMP END PARALLEL DO

       end if
       !
       ! MPI exchange
       !
       call PAR_INTERFACE_NODE_EXCHANGE(ndof2,yy,'SUM','IN MY ZONE','SYNCHRONOUS')

    end if

  end subroutine solver_CSR_smvp

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Matrix vector product
  !> @details Multiply a non symmetric matrix stored in BCSR by a vector
  !>          YY = A XX
  !>          INPUT
  !>          NBNODES .... Number of equations
  !>            NBVAR ...... Number of variables
  !>            AN ......... Matrix
  !>            JA ......... List of elements
  !>            IA ......... Pointer to list of elements
  !>            XX ......... Vector
  !>          OUTPUT
  !>            YY ......... result vector
  !>
  !----------------------------------------------------------------------

  subroutine solver_COO_smvp(nbnodes,nnz,ndof1,ndof2,an,rows,cols,xx,yy)
    integer(ip),  intent(in)           :: nbnodes
    integer(ip),  intent(in)           :: nnz
    integer(ip),  intent(in)           :: ndof1
    integer(ip),  intent(in)           :: ndof2
    real(rp),     intent(in)           :: an(ndof1,ndof2,*)
    integer(ip),  intent(in)           :: rows(*),cols(*)
    real(rp),     intent(in)           :: xx(ndof1,*)
    real(rp),     intent(out)          :: yy(ndof2,*)
    integer(ip)                        :: ii,jj,kk,ll,col,iz
    real(rp)                           :: raux,raux1,raux2,raux3

    if( INOTMASTER ) then
       !
       ! Initialization
       !
       do ii = 1,nbnodes
          yy(1:ndof2,ii) = 0.0_rp
       end do

       if( ndof1 == 1 .and. ndof2 == 1 ) then
          !
          ! NBVAR=1
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)          &
          !$OMP DEFAULT  ( NONE )                      &
          !$OMP PRIVATE  ( ii, jj, iz )                &
          !$OMP SHARED   ( xx, yy, an, nnz, rows, cols )

          do iz = 1,nnz
             ii       = rows(iz)
             jj       = cols(iz)
             yy(1,ii) = yy(1,ii) + an(1,1,iz) * xx(1,jj)
          end do
          !$OMP END PARALLEL DO

       else if( ndof1 == 2 .and. ndof2 == 2 ) then
          !
          ! NBVAR=2
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)            &
          !$OMP DEFAULT  ( NONE )                        &
          !$OMP PRIVATE  ( ii, jj, iz, raux1, raux2    ) &
          !$OMP SHARED   ( xx, yy, an, nnz, rows, cols )
          do iz = 1,nnz
             ii       = rows(iz)
             jj       = cols(iz)
             raux1    = xx(1,jj)
             raux2    = xx(2,jj)
             yy(1,ii) = yy(1,ii) + an(1,1,iz) * raux1
             yy(1,ii) = yy(1,ii) + an(2,1,iz) * raux2
             yy(2,ii) = yy(2,ii) + an(1,2,iz) * raux1
             yy(2,ii) = yy(2,ii) + an(2,2,iz) * raux2
          end do
          !$OMP END PARALLEL DO

       else if( ndof1 == 3 .and. ndof2 == 3 ) then
          !
          ! NBVAR=3
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)                &
          !$OMP DEFAULT  ( NONE )                            &
          !$OMP PRIVATE  ( ii, jj, iz, raux1, raux2, raux3 ) &
          !$OMP SHARED   ( xx, yy, an, nnz, rows, cols     )
          do iz = 1,nnz
             ii       = rows(iz)
             jj       = cols(iz)
             raux1    = xx(1,jj)
             raux2    = xx(2,jj)
             raux3    = xx(3,jj)
             yy(1,ii) = yy(1,ii) + an(1,1,iz) * raux1
             yy(1,ii) = yy(1,ii) + an(2,1,iz) * raux2
             yy(1,ii) = yy(1,ii) + an(3,1,iz) * raux3
             yy(2,ii) = yy(2,ii) + an(1,2,iz) * raux1
             yy(2,ii) = yy(2,ii) + an(2,2,iz) * raux2
             yy(2,ii) = yy(2,ii) + an(3,2,iz) * raux3
             yy(3,ii) = yy(3,ii) + an(1,3,iz) * raux1
             yy(3,ii) = yy(3,ii) + an(2,3,iz) * raux2
             yy(3,ii) = yy(3,ii) + an(3,3,iz) * raux3
          end do
          !$OMP END PARALLEL DO

       else if( ndof1 == 4 .and. ndof2 == 4 ) then
          !
          ! NBVAR=4
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)            &
          !$OMP DEFAULT  ( NONE )                        &
          !$OMP PRIVATE  ( ii, jj, iz, raux            ) &
          !$OMP SHARED   ( xx, yy, an, nnz, rows, cols )
          do iz = 1,nnz
             ii       = rows(iz)
             jj       = cols(iz)
             raux     = xx(1,jj)
             yy(1,ii) = yy(1,ii) + an(1,1,iz) * raux
             yy(2,ii) = yy(2,ii) + an(1,2,iz) * raux
             yy(3,ii) = yy(3,ii) + an(1,3,iz) * raux
             yy(4,ii) = yy(4,ii) + an(1,4,iz) * raux
             raux     = xx(2,jj)
             yy(1,ii) = yy(1,ii) + an(2,1,iz) * raux
             yy(2,ii) = yy(2,ii) + an(2,2,iz) * raux
             yy(3,ii) = yy(3,ii) + an(2,3,iz) * raux
             yy(4,ii) = yy(4,ii) + an(2,4,iz) * raux
             raux     = xx(3,jj)
             yy(1,ii) = yy(1,ii) + an(3,1,iz) * raux
             yy(2,ii) = yy(2,ii) + an(3,2,iz) * raux
             yy(3,ii) = yy(3,ii) + an(3,3,iz) * raux
             yy(4,ii) = yy(4,ii) + an(3,4,iz) * raux
             raux     = xx(4,jj)
             yy(1,ii) = yy(1,ii) + an(4,1,iz) * raux
             yy(2,ii) = yy(2,ii) + an(4,2,iz) * raux
             yy(3,ii) = yy(3,ii) + an(4,3,iz) * raux
             yy(4,ii) = yy(4,ii) + an(4,4,iz) * raux
          end do
          !$OMP END PARALLEL DO

       else
          !
          ! NBVAR = whatever
          !
          !$OMP PARALLEL DO SCHEDULE (STATIC)                        &
          !$OMP DEFAULT  ( NONE )                                    &
          !$OMP PRIVATE  ( ii, jj, ll, iz )                          &
          !$OMP SHARED   ( xx, yy, an, ndof1, ndof2, nnz, rows, cols )
          do iz = 1,nnz
             ii = rows(iz)
             jj = cols(iz)
             do ll = 1,ndof1
                yy(1:ndof2,ii) = yy(1:ndof2,ii) + an(ll,1:ndof2,iz) * xx(ll,jj)
             end do
          end do
          !$OMP END PARALLEL DO

       end if
       !
       ! MPI exchange
       !
       call PAR_INTERFACE_NODE_EXCHANGE(ndof2,yy,'SUM','IN MY ZONE','SYNCHRONOUS')

    end if

  end subroutine solver_COO_smvp

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Solve a system including pre and postprocess
  !> @details Solve a system:
  !>          1. Preprocess: periodicity, reaction force, Neumann b.c.
  !>             and Dirichlet b.c. on A,b
  !>          2. Solve Au=b
  !>          3. Postprocess: reaction
  !>
  !----------------------------------------------------------------------

  subroutine solver_solve(solve,A,b,u,M)
    type(soltyp),      intent(inout), pointer  :: solve(:)
    real(rp),          intent(inout)           :: A(*)
    real(rp),          intent(inout)           :: b(*)
    real(rp),          intent(inout)           :: u(*)
    real(rp),          intent(inout), optional :: M(*)
    
    if( solve(1) % kfl_version == 1 ) call solver_preprocess(solve,A,b,u)
    if( present(M) ) then
       call solver_solve_system(solve,A,b,u,M)
    else
       call solver_solve_system(solve,A,b,u)
    end if
    if( solve(1) % kfl_version == 1 ) call solver_postprocess(solve,A,b,u)

  end subroutine solver_solve

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Solve a system
  !> @details Solve the system Au=b. M can be a preconditioner
  !>
  !----------------------------------------------------------------------

  subroutine solver_solve_system(solve,A,b,u,M)
    type(soltyp),      intent(in),    pointer  :: solve(:)
    real(rp),          intent(inout)           :: A(*)
    real(rp),          intent(inout)           :: b(*)
    real(rp),          intent(inout)           :: u(*)
    real(rp),          intent(inout), optional :: M(*)
    real(rp)                                   :: dummr(2)

    if( present(M) ) then
       call solver(b,u,A,M)
    else
       call solver(b,u,A,dummr)
    end if

  end subroutine solver_solve_system

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Solver preprocess
  !> @details Preprocess the matrix and RHS A,b to account for:\n
  !>          1. Impose periodicity \n
  !>          2. Save matrix and RHS for reaction force \n
  !>          3. Impose Neumann boundary conditionsAssemble reaction \n
  !>          4. Impose Dirichlet boundary conditions \n
  !>
  !----------------------------------------------------------------------

  subroutine solver_preprocess(solve,A,b,u,P)
    
    type(soltyp), intent(inout),    pointer  :: solve(:)
    real(rp),     intent(inout)              :: A(*)
    real(rp),     intent(inout)              :: b(*)
    real(rp),     intent(inout)              :: u(*)
    real(rp),     intent(inout),    optional :: P(*)
    integer(ip)                              :: ndof1,ndof2
    integer(ip)                              :: pointer_a11
    integer(ip)                              :: pointer_a12
    integer(ip)                              :: pointer_a21
    integer(ip)                              :: pointer_a22
    integer(ip)                              :: pointer_b1
    integer(ip)                              :: pointer_b2
    integer(ip)                              :: num_blocks
    real(rp)                                 :: dummr(3)
    character(14)                            :: matrix_rhs
    
    if( INOTMASTER ) then

       num_blocks = solve(1) % num_blocks
       if( solve(1) % kfl_algso == SOL_SOLVER_RICHARDSON ) then
          matrix_rhs = 'RHS'
       else
          matrix_rhs = 'MATRIX AND RHS'
       end if

       !-----------------------------------------------------------------
       !
       ! 1. Impose periodicity
       ! 2. Save matrix and RHS for reaction force
       ! 3. Assemble reaction
       ! 4. Impose boundary conditions
       !
       !-----------------------------------------------------------------

       if( num_blocks == 1 ) then

          ndof1 = solve(1) % ndofn
          !
          ! 1. Periodicity
          !
          call solver_periodicity(trim(matrix_rhs),solve,ndof1,ndof1,A,b,u)
          !
          ! 2. Reaction force: save matrix and RHS
          !
          call solver_reaction_force_1by1('SAVE SYSTEM',solve,ndof1,A,b,u)
          !
          ! 3. Assemble reaction force
          !
          call solver_reaction_force_1by1('ASSEMBLE REACTION',solve,ndof1,A,b,u)
          !
          ! 4. Impose Dirichlet
          !
          !allocate(solver_invdiag(npoin))
          !call solver_diagonal(solve,A,solver_invdiag)

          call solver_preprocess_dirichlet_1by1(trim(matrix_rhs),solve,ndof1,A,b,u)

       else if( num_blocks == 2 ) then
          !
          ! Check dimensions
          !
          ndof1       = solve(1) % block_dimensions(1)
          ndof2       = solve(1) % block_dimensions(2)
          pointer_a11 = 1
          pointer_a12 = pointer_a11 + nzdom * ndof1 * ndof1
          pointer_a21 = pointer_a12 + nzdom * ndof1 * ndof2
          pointer_a22 = pointer_a21 + nzdom * ndof2 * ndof1
          pointer_b1  = 1
          pointer_b2  = pointer_b1  + npoin * ndof1
          !
          ! 1. Periodicity
          !
          call solver_periodicity('MATRIX AND RHS',solve,ndof1,ndof1,A(pointer_a11),b(pointer_b1),u(pointer_b1))
          call solver_periodicity('MATRIX',        solve,ndof1,ndof2,A(pointer_a12),b(pointer_b1),u(pointer_b1))
          call solver_periodicity('MATRIX',        solve,ndof2,ndof1,A(pointer_a21),b(pointer_b2),u(pointer_b2))
          call solver_periodicity('MATRIX AND RHS',solve,ndof2,ndof2,A(pointer_a22),b(pointer_b2),u(pointer_b2))
          if( present(P) ) then
             call solver_periodicity('MATRIX',solve,ndof2,ndof2,P,dummr,dummr,'KEEP DIAGONAL')
          end if
          !
          ! 2. Reaction force: save matrix and RHS
          !
          call solver_reaction_force_2by2(&
               & 'SAVE SYSTEM',solve,ndof1,ndof2,&
               & A(pointer_a11),A(pointer_a12),A(pointer_a21),A(pointer_a22),&
               & b(pointer_b1), b(pointer_b2), u(pointer_b1), u(pointer_b2))
          !
          ! 3. Assemble reaction force
          !
          call solver_reaction_force_2by2(&
               & 'ASSEMBLE REACTION',solve,ndof1,ndof2,&
               & A(pointer_a11),A(pointer_a12),A(pointer_a21),A(pointer_a22),&
               & b(pointer_b1), b(pointer_b2), u(pointer_b1), u(pointer_b2))
          !
          ! 4. Dirichlet
          !
          ! Not coded for now: When periodicity and local axes, take basis of master for slaves
          !
       else

          call runend('SOLVER_PREPROCESS: NOT READY FOR THIS KIND OF BLOCK MATRIX')

       end if
       
    end if
    
  end subroutine solver_preprocess

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Initialization
  !> @details Initialize matrix and RHS
  !>
  !----------------------------------------------------------------------

  subroutine solver_initialize_matrix_and_rhs_1(solve,A,b,P)
    type(soltyp), intent(inout),    pointer           :: solve(:)
    real(rp),     intent(inout),    pointer           :: A(:)
    real(rp),     intent(inout),    pointer, optional :: b(:)
    real(rp),     intent(inout),    pointer, optional :: P(:)
    integer(ip)                                       :: ndof1,ndof2
    integer(ip)                                       :: num_blocks
    integer(ip)                                       :: nzd11,nzd12,nzd21,nzd22
    integer(ip)                                       :: nrhs1,nrhs2,nrhsa
    integer(ip)                                       :: nzdal,ii

    if( INOTMASTER ) then

       num_blocks = solve(1) % num_blocks

       if( num_blocks == 1 ) then

          !--------------------------------------------------------------
          !
          ! 1 block
          !
          !--------------------------------------------------------------

          if( solve(1) % nzmat > size(A,KIND=ip) ) call runend('SOLVER_INITIALIZE_MATRIX_AND_RHS: WRONG MATRIX DIMENSION')
          do ii = 1,solve(1) % nzmat
             A(ii) = 0.0_rp
          end do
          
          if( present(b) ) then
             if( solve(1) % nzrhs > size(b,KIND=ip) ) call runend('SOLVER_INITIALIZE_MATRIX_AND_RHS: WRONG RHS DIMENSION')
             do ii = 1,solve(1) % nzrhs
                b(ii) = 0.0_rp
             end do
          end if

       else if( num_blocks == 2 ) then

          !--------------------------------------------------------------
          !
          ! 2 blocks
          !
          !--------------------------------------------------------------

          ndof1 = solve(1) % block_dimensions(1)
          ndof2 = solve(1) % block_dimensions(2)
          !
          ! Non-zero coefficients of each block
          !
          nzd11 = solve(1) % nzmat
          nzd12 = solve(1) % nzmat / ndof1 * ndof2
          nzd21 = nzd12
          nzd22 = solve(2) % nzmat
          nzdal = nzd11 + nzd12 + nzd21 + nzd22 ! Non-zero coefficients of commplete matrix
          !
          ! Size of unknowns of each block
          !
          nrhs1 = solve(1) % nzrhs
          nrhs2 = solve(2) % nzrhs
          nrhsa = nrhs1 + nrhs2
          !
          ! Check errors
          !
          if( nzdal > size(A,KIND=ip) ) call runend('SOLVER_INITIALIZE_MATRIX_AND_RHS: WRONG MATRIX DIMENSION')
          if( nrhsa > size(b,KIND=ip) ) call runend('SOLVER_INITIALIZE_MATRIX_AND_RHS: WRONG RHS DIMENSION')
          !
          ! Initialize system
          !
          do ii = 1,nzdal
             A(ii) = 0.0_rp
          end do
          
          if( present(b) ) then
             do ii = 1,nrhsa
                b(ii) = 0.0_rp
             end do
          end if
          
          if( present(P) ) then
             if( solve(2) % nzmat > size(P,KIND=ip) ) call runend('SOLVER_INITIALIZE_MATRIX_AND_RHS: WRONG MATRIX P DIMENSION')
             do ii = 1,solve(2) % nzmat
                P(ii) = 0.0_rp
             end do
          end if

       else

          call runend('SOLVER_INITIALIZE_MATRIX_AND_RHS: NOT READY FOR THIS KIND OF BLOCK MATRIX')

       end if

    end if

  end subroutine solver_initialize_matrix_and_rhs_1

    subroutine solver_initialize_matrix_and_rhs_S(solve,A,b,P)
    type(soltyp), intent(inout),    pointer           :: solve
    real(rp),     intent(inout),    pointer           :: A(:)
    real(rp),     intent(inout),    pointer, optional :: b(:)
    real(rp),     intent(inout),    pointer, optional :: P(:)
    integer(ip)                                       :: ii

    if( INOTMASTER ) then

       if( solve % num_blocks == 1 ) then

          !--------------------------------------------------------------
          !
          ! 1 block
          !
          !--------------------------------------------------------------

          if( solve % nzmat > size(A,KIND=ip) ) call runend('SOLVER_INITIALIZE_MATRIX_AND_RHS: WRONG MATRIX DIMENSION')
          do ii = 1,solve % nzmat
             A(ii) = 0.0_rp
          end do

          if( present(b) ) then
             if( solve % nzrhs > size(b,KIND=ip) ) call runend('SOLVER_INITIALIZE_MATRIX_AND_RHS: WRONG RHS DIMENSION')
             do ii = 1,solve % nzrhs
                b(ii) = 0.0_rp
             end do
          end if
       else
          
          call runend('SOLVER_INITIALIZE_MATRIX_AND_RHS: NOT READY FOR THIS KIND OF BLOCK MATRIX')

       end if

    end if

  end subroutine solver_initialize_matrix_and_rhs_S

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    07/10/2014
  !> @brief   Solver postprocess
  !> @details Posprocess the matrix RHS and u to account for: \n
  !>          1. Periodicity: recover solution on slave nodes \n
  !>          2. Compute reaction \n
  !>
  !----------------------------------------------------------------------

  subroutine solver_postprocess(solve,A,b,u)

    type(soltyp),      intent(inout),    pointer  :: solve(:)
    real(rp),          intent(inout)           :: A(*)
    real(rp),          intent(inout)           :: b(*)
    real(rp),          intent(inout)           :: u(*)
    integer(ip)                                :: ndof1,ndof2
    integer(ip)                                :: pointer_a11
    integer(ip)                                :: pointer_a12
    integer(ip)                                :: pointer_a21
    integer(ip)                                :: pointer_a22
    integer(ip)                                :: pointer_b1
    integer(ip)                                :: pointer_b2
    integer(ip)                                :: num_blocks
    real(rp)                                   :: dummr(3)
    character(14)                              :: matrix_rhs

    if( INOTMASTER ) then

       num_blocks = solve(1) % num_blocks
       if( solve(1) % kfl_algso == SOL_SOLVER_RICHARDSON ) then
          matrix_rhs = 'RHS'
       else
          matrix_rhs = 'MATRIX AND RHS'
       end if

       !-----------------------------------------------------------------
       !
       ! 1. Periodicity: recover solution on slave nodes
       ! 2. Compute reaction
       !
       !-----------------------------------------------------------------

       if( num_blocks == 1 ) then

          ndof1 = solve(1) % ndofn
          !
          ! 1. Periodicity
          !
          call solver_periodicity('SOLUTION',solve,ndof1,ndof1,A,b,u)
          !
          ! 2. Compute reaction force
          !
          call solver_reaction_force_1by1('REACTION FORCE',solve,ndof1,A,b,u)

       else if( num_blocks == 2 ) then
          !
          ! Define block dimensions and pointers
          !
          ndof1       = solve(1) % block_dimensions(1)
          ndof2       = solve(1) % block_dimensions(2)
          pointer_a11 = 1
          pointer_a12 = pointer_a11 + nzdom * ndof1 * ndof1
          pointer_a21 = pointer_a12 + nzdom * ndof1 * ndof2
          pointer_a22 = pointer_a21 + nzdom * ndof2 * ndof1
          pointer_b1  = 1
          pointer_b2  = pointer_b1  + npoin * ndof1
          !
          ! 1. Periodicity
          !
          call solver_periodicity('SOLUTION',solve,ndof1,ndof1,dummr,dummr,u(pointer_b1))
          call solver_periodicity('SOLUTION',solve,ndof2,ndof2,dummr,dummr,u(pointer_b2))
          !
          ! 2. Compute reaction force
          !
          call solver_reaction_force_2by2(&
               & 'REACTION FORCE',solve,ndof1,ndof2,&
               & A(pointer_a11),A(pointer_a12),A(pointer_a21),A(pointer_a22),&
               & b(pointer_b1), b(pointer_b2), u(pointer_b1), u(pointer_b2))

       else

          call runend('SOLVER_POSTPROCESS: NOT READY FOR THIS KIND OF BLOCK MATRIX')

       end if

    end if

  end subroutine solver_postprocess

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    03/10/2014
  !> @brief   Impose periodicity
  !> @details Periodicity is imposed in two steps, as explaine below
  !>          \verbatim
  !>
  !>               1 M 3 4 S 6 7 8             1 M 3 4 S 6 7 8             1 M 3 4 S 6 7 8
  !>            +-                -+        +-                -+        +-                -+
  !>          1 |  x x             |      1 |  x x             |      1 |  x x             |
  !>          M |  x x x   ||      |      M |  x x x           |      M |  x x x X X X <=  |
  !>          3 |    x x x \/      |      3 |    x x x         |      3 |    x x x         |
  !>          4 |      x x x       |  =>  4 | => X x x 0       |  =>  4 |      x x x       |
  !>          S |        x x x     |      S |        x x x     |      S |        0 0 0 <=  |
  !>          6 |     <--- x x x   |      6 | => X     0 x x   |      6 |          x x x   |
  !>          7 |     move   x x x |      7 |            x x x |      7 |            x x x |
  !>          8 |              x x |      8 |              x x |      8 |              x x |
  !>            +-                -+        +-                -+        +-                -+
  !>
  !>            Original matrix             For nodes who are not       Move slave's rows to
  !>                                        slaves, move slave's        their master's rows
  !>                                        column to master's one
  !>           \endverbatim
  !>
  !>           Comments:
  !>           - The master is always in the subdomains where it's slave is
  !>           - The master is in the inter-partition boundary list, slave may not
  !>
  !----------------------------------------------------------------------

  subroutine solver_periodicity(wtask,solve,ndof1,ndof2,A,b,u,message)
    character(*),      intent(in)             :: wtask
    !type(soltyp),      intent(in),    pointer :: solve(:)
    type(soltyp),      intent(in)             :: solve(*)
    integer(ip),       intent(in)             :: ndof1
    integer(ip),       intent(in)             :: ndof2
    real(rp),          intent(inout)          :: A(ndof2,ndof1,nzdom)
    real(rp),          intent(inout),optional :: b(ndof1,npoin)
    real(rp),          intent(inout),optional :: u(ndof1,npoin)
    character(*),      intent(in),   optional :: message
    integer(ip)                               :: izdom,jzdom,ifoun
    integer(ip)                               :: ipoin,jpoin,kpoin,jj,ii
    integer(ip)                               :: iperi,kzdom,lpoin
    integer(ip),                      pointer :: list_of_masters(:)
    logical(lg),                      pointer :: in_my_zone(:)
    logical(lg)                               :: periodicity_matrix
    logical(lg)                               :: periodicity_rhs
    logical(lg)                               :: periodicity_solution

    if( IMASTER )    return
    if( nperi == 0 ) return
    nullify( list_of_masters )
    nullify( in_my_zone )
    !
    ! Decide what to do
    !
    select case ( trim(wtask) )
    case ( 'RHS' )
       periodicity_matrix   = .false.
       periodicity_rhs      = .true.
       periodicity_solution = .false.
    case ( 'MATRIX' )
       periodicity_matrix   = .true.
       periodicity_rhs      = .false.
       periodicity_solution = .false.
    case ( 'SOLUTION' )
       periodicity_matrix   = .false.
       periodicity_rhs      = .false.
       periodicity_solution = .true.
    case ( 'MATRIX AND RHS' )
       periodicity_matrix   = .true.
       periodicity_rhs      = .true.
       periodicity_solution = .false.
    case default
       call runend('PERIODICITY: DO NOT KNOW WHAT TO DO')
    end select

    if( periodicity_rhs .and. ( .not. periodicity_matrix .and. .not. periodicity_solution ) .and. present(b) ) then

       !--------------------------------------------------------------
       !
       ! Only RHS
       !
       !--------------------------------------------------------------

       if( nzone > 1 ) then
          !do kpoin = 1,npoiz(current_zone)
          !   jpoin = lpoiz(current_zone) % l(kpoin)
          !   ipoin = lmast(jpoin)
          !   if( ipoin > 0 ) then
          !      if( in_my_zone(ipoin) .and. in_my_zone(jpoin) ) then
          !         b(1:ndof1,ipoin) = b(1:ndof1,ipoin) + b(1:ndof1,jpoin)
          !         b(1:ndof1,jpoin) = b(1:ndof1,ipoin)
          !      end if
          !   end if
          !end do
          call memory_alloca(memit,'IN_MY_ZONE','solver_periodicity',in_my_zone,npoin)
          do ipoin = 1,npoin
             in_my_zone(ipoin) = .true.
          end do
          do iperi = 1,nperi
             ipoin = lperi(1,iperi)
             jpoin = lperi(2,iperi)
             if ( ipoin > 0 .and. jpoin > 0 ) then
                if( in_my_zone(ipoin) .and. in_my_zone(jpoin) ) then
                   b(1:ndof1,ipoin) = b(1:ndof1,ipoin) + b(1:ndof1,jpoin)
                   b(1:ndof1,jpoin) = b(1:ndof1,ipoin)
                end if
             end if
          end do
          call memory_deallo(memit,'IN_MY_ZONE','solver_periodicity',in_my_zone)
       else
          !do jpoin = 1,npoin
          !   ipoin = lmast(jpoin)
          !   if( ipoin > 0 ) then
          !      b(1:ndof1,ipoin) = b(1:ndof1,ipoin) + b(1:ndof1,jpoin)
          !      b(1:ndof1,jpoin) = b(1:ndof1,ipoin)
          !   end if
          !end do
          do iperi = 1,nperi
             ipoin = lperi(1,iperi)
             jpoin = lperi(2,iperi)
             if( ipoin > 0 .and. jpoin > 0 ) then
                b(1:ndof1,ipoin) = b(1:ndof1,ipoin) + b(1:ndof1,jpoin)
                b(1:ndof1,jpoin) = b(1:ndof1,ipoin)
             end if
          end do
       end if

    else if( periodicity_matrix ) then

       !--------------------------------------------------------------
       !
       ! Matrix and/or RHS
       !
       !--------------------------------------------------------------
       !
       ! Check errors
       !
       if( solve(1) % kfl_algso == 0 ) then
          call runend('PERIODICITY ONLY POSSIBLE WITH SPARSE SOLVERS')
       else if( solve(1) % kfl_symme == 1 ) then
          call runend('MATRIX PERIODICITY ONLY POSSIBLE WITH ASYMMETRIC SOLVERS')
       end if
       !
       ! Put all slaves' columns into master's column
       ! KPOIN is neither master nor slave
       ! LPERI(2,:) = JPOIN is slave
       ! LPERI(1,:) = IPOIN is JPOIN's master
       ! The master is in all the partitions where slave is present, so the entire
       ! row of the slave can be assembled into the master's row
       !
       call memory_alloca(memit,'IN_MY_ZONE',     'solver_periodicity',in_my_zone,npoin)
       call memory_alloca(memit,'LIST_OF_MASTERS','solver_periodicity',list_of_masters,npoin)
       do ipoin = 1,npoin
          in_my_zone(ipoin) = .true.
       end do
       do iperi = 1,nperi
          ipoin = lperi(1,iperi)
          jpoin = lperi(2,iperi)
          if ( ipoin <= 0 .or. jpoin <= 0 ) then
             in_my_zone(abs(ipoin)) = .false.
          end if
       end do

       do iperi = 1,nperi
          ipoin = lperi(1,iperi)
          jpoin = lperi(2,iperi)
          if ( ipoin > 0 .and. jpoin > 0 ) then
             if( in_my_zone(ipoin) .and. in_my_zone(jpoin) ) then
                list_of_masters(jpoin) =  ipoin
                !list_of_masters(ipoin) = -1
             end if
             !list_of_masters(jpoin) = ipoin ! Ask eva
          end if
       end do
       !
       ! Move slave's column to master's column
       !
       do kpoin = 1,npoin

          if( list_of_masters(kpoin) == 0 ) then
             do kzdom = r_dom(kpoin),r_dom(kpoin+1)-1
                jpoin = c_dom(kzdom)
                if( list_of_masters(jpoin) > 0 .and. in_my_zone(jpoin) ) then
                   ipoin = list_of_masters(jpoin)
                   ifoun = 0
                   izdom2: do izdom = r_dom(kpoin),r_dom(kpoin+1)-1
                      if( c_dom(izdom) == ipoin ) then
                         ifoun = 1

                         do ii = 1,ndof2
                            do jj = 1,ndof1
                               A(ii,jj,izdom) = A(ii,jj,izdom) + A(ii,jj,kzdom)
                               A(ii,jj,kzdom) = 0.0_rp
                            end do
                         end do

                         exit izdom2
                      end if
                   end do izdom2
                   if( ifoun == 0 ) then
                      print*,'NOT FOUND NODE, SLAVE, MASTER A=',kpoin,jpoin,ipoin
                      call runend('PRESOL 2: NODE NOT FOUND')
                   end if
                end if

             end do

          end if
       end do
       !
       ! Put all slaves' rows into master's row
       ! Put slave row to zero
       ! JPOIN:  slave
       ! IPOIN:  master
       ! Slave:  coef. JZDOM for JPOIN-KPOIN
       ! Master: coef. IZDOM for IPOIN-KPOIN
       !
       do jpoin = 1,npoin
          ipoin = list_of_masters(jpoin)
          if( ipoin > 0 ) then
             if( periodicity_rhs .and. present(b) ) then
                b(1:ndof1,ipoin) = b(1:ndof1,ipoin) + b(1:ndof1,jpoin)
                b(1:ndof1,jpoin) = 0.0_rp
             end if
             do jzdom = r_dom(jpoin),r_dom(jpoin+1)-1
                kpoin = c_dom(jzdom)
                if( in_my_zone(kpoin) ) then
                   if( list_of_masters(kpoin) > 0 ) kpoin = list_of_masters(kpoin)
                   ifoun = 0
                   izdom1: do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                      if( c_dom(izdom) == kpoin ) then
                         ifoun = 1
                         do ii = 1,ndof2
                            do jj = 1,ndof1
                               A(ii,jj,izdom) = A(ii,jj,izdom) + A(ii,jj,jzdom)
                               A(ii,jj,jzdom) = 0.0_rp
                            end do
                         end do
                         exit izdom1
                      end if
                   end do izdom1
                   if( ifoun == 0 .and. ipoin /= kpoin ) then
                      !print*,'NOT FOUND NODE, SLAVE, MASTER B=',kfl_paral,lninv_loc(kpoin),lninv_loc(jpoin),lninv_loc(ipoin)
                      call runend('PRESOL 1: NODE NOT FOUND')
                   end if
                end if
             end do
          end if
       end do

       if( present(message) ) then
          if( trim(message) == 'KEEP DIAGONAL') then
             do iperi = 1,nperi
                ipoin = lperi(2,iperi)
                jpoin = lperi(1,iperi)
                if ( ipoin > 0 .and. jpoin > 0 ) then
                   loop_slave: do kzdom = r_dom(ipoin),r_dom(ipoin+1)-1
                      if( c_dom(kzdom) == ipoin ) then
                         do ii = 1,ndof1
                            A(ii,ii,kzdom) = 1.0_rp
                         end do
                         exit loop_slave
                      end if
                   end do loop_slave
                end if
             end do
          end if
       end if
       !do jpoin = 1,npoin
       !   if( list_of_masters(jpoin) > 0 ) then
       !      do jzdom = r_dom(jpoin),r_dom(jpoin+1)-1
       !         kpoin = c_dom(jzdom)
       !         if( kpoin == jpoin ) then
       !            do ii = 1,ndof2
       !               do jj = 1,ndof1
       !                  if( ii /= jj ) then
       !                     A(ii,jj,jzdom) = 0.0_rp
       !                  end if
       !               end do
       !            end do
       !         else
       !            A(1:ndof2,1:ndof1,jzdom) = 0.0_rp
       !         end if
       !      end do
       !   end if
       !end do

       call memory_deallo(memit,'IN_MY_ZONE',     'solver_periodicity',in_my_zone)
       call memory_deallo(memit,'LIST_OF_MASTERS','solver_periodicity',list_of_masters)

    else if( periodicity_solution .and. present(u) ) then

       !--------------------------------------------------------------
       !
       ! Solution
       !
       !--------------------------------------------------------------

       if( nzone > 1 ) then
          call memory_alloca(memit,'IN_MY_ZONE','solver_periodicity',in_my_zone,npoin)
          do ipoin = 1,npoin
             in_my_zone(ipoin) = .true.
          end do
          do iperi = 1,nperi
             jpoin = lperi(2,iperi)
             ipoin = lperi(1,iperi)
             if( ipoin > 0 .and. jpoin > 0 ) then
                if( in_my_zone(ipoin) .and. in_my_zone(jpoin) ) then
                   u(1:ndof1,jpoin) = u(1:ndof1,ipoin)
                end if
             end if
          end do
          call memory_deallo(memit,'IN_MY_ZONE','solver_periodicity',in_my_zone)
       else          
          do iperi = 1,nperi
             jpoin = lperi(2,iperi)
             ipoin = lperi(1,iperi)
             if( ipoin > 0 .and. jpoin > 0 ) u(1:ndof1,jpoin) = u(1:ndof1,ipoin)
          end do
       end if

    end if

  end subroutine solver_periodicity

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    03/10/2014
  !> @brief   Impose perdiodicity on RHS
  !> @details Impose periodicity on the RHS of a system:
  !>
  !>          B(MASTER) = B(MASTER) + B(SLAVE)
  !>          B(SLAVE)  = B(MASTER)
  !>
  !----------------------------------------------------------------------

  subroutine solver_periodicity_rhs(ndof1,b)
    integer(ip),       intent(in)    :: ndof1          !< # dof
    real(rp),          intent(inout) :: b(ndof1,npoin) !< RHS
    integer(ip)                      :: ipoin,jpoin
    integer(ip)                      :: iperi

    if( INOTMASTER ) then

       do iperi = 1,nperi
          ipoin = lperi(1,iperi) ! Master
          jpoin = lperi(2,iperi) ! Slave
          if( ipoin > 0 .and. jpoin > 0 ) then
             b(1:ndof1,ipoin) = b(1:ndof1,ipoin) + b(1:ndof1,jpoin)
             b(1:ndof1,jpoin) = b(1:ndof1,ipoin)
          end if
       end do
    end if

  end subroutine solver_periodicity_rhs

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    03/10/2014
  !> @brief   Impose perdiodicity on solution
  !> @details Impose periodicity on the solution U of a periodic system:
  !>
  !>          U(SLAVE)  = U(MASTER)
  !>
  !----------------------------------------------------------------------

  subroutine solver_periodicity_solution(ndof1,u)
    integer(ip),       intent(in)    :: ndof1          !< # dof
    real(rp),          intent(inout) :: u(ndof1,npoin) !< Solution
    integer(ip)                      :: ipoin,jpoin
    integer(ip)                      :: iperi

    if( INOTMASTER ) then

       do iperi = 1,nperi
          ipoin = lperi(1,iperi)    ! Master
          jpoin = lperi(2,iperi)    ! Slave
          if( ipoin > 0 .and. jpoin > 0 ) u(1:ndof1,jpoin) = u(1:ndof1,ipoin)
       end do

    end if

  end subroutine solver_periodicity_solution

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    03/10/2014
  !> @brief   Save matrix and RHS to compute reaction forces
  !> @details Save matrix and RHS to compute reaction forces
  !>
  !----------------------------------------------------------------------

  subroutine solver_reaction_force_1by1(wtask,solve,ndof1,A11,b1,u1)
    use def_domain, only : coord
    character(*),      intent(in)             :: wtask
    type(soltyp),      intent(inout),    pointer :: solve(:)
    integer(ip),       intent(in)             :: ndof1
    real(rp),          intent(inout)          :: A11(ndof1,ndof1,*)
    real(rp),          intent(inout)          :: b1(ndof1,*)
    real(rp),          intent(inout)          :: u1(ndof1,*)
    integer(ip)                               :: izdom,jzdom,lpoin,ipoin
    integer(ip)                               :: idofn,jdofn,jpoin

    if( associated(solve(1) % lpoin_reaction) ) then

       if( solve(1) % num_blocks /= 1 ) call runend('1BY1: WRONG NUMBER OF BLOCKS')

       select case ( trim(wtask) )

       case ( 'SAVE SYSTEM' )
          !
          ! Save system: matrix and RHS at reaction nodes
          !
          do ipoin = 1,npoin
             if( solve(1) % lpoin_reaction(ipoin) ) then
                solve(1) % lpoin_block(ipoin) % block1_num(1) % rhs(1:ndof1) = b1(1:ndof1,ipoin)
                jzdom = 0
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jzdom = jzdom + 1
                   solve(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix(1:ndof1,1:ndof1,jzdom) = A11(1:ndof1,1:ndof1,izdom)
                end do
             end if
          end do

       case ( 'REACTION FORCE' )
          !
          ! Compute reaction force
          !
          if( solve(1) % kfl_react /= 0 ) then

             do ipoin = 1,npoin
                solve(1) % reaction(1:ndof1,ipoin) = 0.0_rp
                if( solve(1) % lpoin_reaction(ipoin) ) then
                   jzdom = 0
                   solve(1) % reaction(1:ndof1,ipoin) = solve(1) % lpoin_block(ipoin) % block1_num(1) % rhs(1:ndof1)
                   do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                      jpoin = c_dom(izdom)
                      jzdom = jzdom + 1
                      do idofn = 1,ndof1
                         do jdofn = 1,ndof1
                            solve(1) % reaction(idofn,ipoin) = solve(1) % reaction(idofn,ipoin) &
                                 & - solve(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix(jdofn,idofn,jzdom) * u1(jdofn,jpoin)
                         end do
                      end do
                   end do
                end if
             end do
             !call pararr('SOL',NPOIN_TYPE,NPOIN,solve(1) % reaction)
             !
             ! Remove this line to see reaction on coupling interefaces
             !
             call PAR_INTERFACE_NODE_EXCHANGE(solve(1) % reaction,'SUM','IN MY ZONE')

          end if

       case ( 'ASSEMBLE REACTION' )
          !
          ! Assemble reaction
          !
          if( solve(1) % kfl_bvnat == 1 ) then
             if( solve(1) % kfl_iffix /= 0 ) then
                do ipoin = 1,npoin
                   do idofn = 1,ndof1
                      if( solve(1) % kfl_fixno(idofn,ipoin) <= 0 ) then
                         b1(idofn,ipoin) = b1(idofn,ipoin) + solve(1) % bvnat(idofn,ipoin)
                      end if
                   end do
                end do
             else
                do ipoin = 1,npoin
                   b1(1:ndof1,ipoin) = b1(1:ndof1,ipoin) + solve(1) % bvnat(1:ndof1,ipoin)
                end do
             end if
          end if

       case default
          !
          ! Do not know what to do
          !
          call runend('1BY1: DOES KNOW WHAT TO DO')

       end select

    end if

  end subroutine solver_reaction_force_1by1

  subroutine solver_reaction_force_2by2(wtask,solve,ndof1,ndof2,A11,A12,A21,A22,b1,b2,u1,u2)
    character(*), intent(in)             :: wtask
    type(soltyp), intent(inout), pointer :: solve(:)
    integer(ip),  intent(in)             :: ndof1
    integer(ip),  intent(in)             :: ndof2
    real(rp),     intent(inout)          :: A11(ndof1,ndof1,*)
    real(rp),     intent(inout)          :: A12(ndof2,ndof1,*)
    real(rp),     intent(inout)          :: A21(ndof1,ndof2,*)
    real(rp),     intent(inout)          :: A22(ndof2,ndof2,*)
    real(rp),     intent(inout)          :: b1(ndof1,*)
    real(rp),     intent(inout)          :: b2(ndof2,*)
    real(rp),     intent(inout)          :: u1(ndof1,*)
    real(rp),     intent(inout)          :: u2(ndof2,*)
    integer(ip)                          :: izdom,jzdom,lpoin,ipoin
    integer(ip)                          :: jpoin,jdofn,idofn

    if( associated(solve(1) % lpoin_reaction) ) then

       if( solve(1) % num_blocks /= 2 ) call runend('2BY2: WRONG NUMBER OF BLOCKS')

       select case ( trim(wtask) )

       case ( 'SAVE SYSTEM' )
          !
          ! Save system: matrix and RHS at reaction nodes
          !
          do ipoin = 1,npoin
             if( solve(1) % lpoin_reaction(ipoin) ) then
                solve(1) % lpoin_block(ipoin) % block1_num(1) % rhs(1:ndof1) = b1(1:ndof1,ipoin)
                solve(1) % lpoin_block(ipoin) % block1_num(2) % rhs(1:ndof2) = b2(1:ndof2,ipoin)
                jzdom = 0
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jzdom = jzdom + 1
                   solve(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix(1:ndof1,1:ndof1,jzdom) = A11(1:ndof1,1:ndof1,izdom)
                   solve(1) % lpoin_block(ipoin) % block2_num(1,2) % matrix(1:ndof2,1:ndof1,jzdom) = A12(1:ndof2,1:ndof1,izdom)
                   solve(1) % lpoin_block(ipoin) % block2_num(2,1) % matrix(1:ndof1,1:ndof2,jzdom) = A21(1:ndof1,1:ndof2,izdom)
                   solve(1) % lpoin_block(ipoin) % block2_num(2,2) % matrix(1:ndof2,1:ndof2,jzdom) = A22(1:ndof2,1:ndof2,izdom)
                end do
             end if
          end do

       case ( 'REACTION FORCE' )
          !
          ! Compute reaction force
          !
          if( solve(1) % kfl_react /= 0 ) then

             do ipoin = 1,npoin
                solve(1) % reaction(1:ndof1,ipoin) = 0.0_rp
                solve(2) % reaction(1:ndof2,ipoin) = 0.0_rp
                if( solve(1) % lpoin_reaction(ipoin) ) then
                   jzdom = 0
                   solve(1) % reaction(1:ndof1,ipoin) = solve(1) % lpoin_block(ipoin) % block1_num(1) % rhs(1:ndof1)
                   solve(2) % reaction(1:ndof2,ipoin) = solve(1) % lpoin_block(ipoin) % block1_num(2) % rhs(1:ndof2)
                   do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                      jpoin = c_dom(izdom)
                      jzdom = jzdom + 1
                      do idofn = 1,ndof1
                         do jdofn = 1,ndof1
                            solve(1) % reaction(idofn,ipoin) = solve(1) % reaction(idofn,ipoin) &
                                 & - solve(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix(jdofn,idofn,jzdom) &
                                 & * u1(jdofn,jpoin)
                         end do
                         do jdofn = 1,ndof2
                            solve(1) % reaction(idofn,ipoin) = solve(1) % reaction(idofn,ipoin) &
                                 & - solve(1) % lpoin_block(ipoin) % block2_num(1,2) % matrix(jdofn,idofn,jzdom) &
                                 & * u2(jdofn,jpoin)
                         end do
                      end do
                      do idofn = 1,ndof2
                         do jdofn = 1,ndof1
                            solve(2) % reaction(idofn,ipoin) = solve(2) % reaction(idofn,ipoin) &
                                 & - solve(1) % lpoin_block(ipoin) % block2_num(2,1) % matrix(jdofn,idofn,jzdom) &
                                 & * u1(jdofn,jpoin)
                         end do
                         do jdofn = 1,ndof2
                            solve(2) % reaction(idofn,ipoin) = solve(2) % reaction(idofn,ipoin) &
                                 & - solve(1) % lpoin_block(ipoin) % block2_num(2,2) % matrix(jdofn,idofn,jzdom) &
                                 & * u2(jdofn,jpoin)
                         end do
                      end do
                   end do
                end if
             end do
             !call pararr('SLX',NPOIN_TYPE,NPOIN,solve(1) % reaction)
             !call pararr('SLX',NPOIN_TYPE,NPOIN,solve(2) % reaction)
             call PAR_INTERFACE_NODE_EXCHANGE(solve(1) % reaction,'SUM','IN MY ZONE')
             call PAR_INTERFACE_NODE_EXCHANGE(solve(2) % reaction,'SUM','IN MY ZONE')

          end if

       case ( 'ASSEMBLE REACTION' )
          !
          ! Assemble reaction
          !
          if( solve(1) % kfl_bvnat == 1 ) then
             if( solve(1) % kfl_iffix /= 0 ) then
                do ipoin = 1,npoin
                   do idofn = 1,ndof1
                      if( solve(1) % block_array(1) % kfl_fixno(idofn,ipoin) <= 0 ) then
                         b1(idofn,ipoin) = b1(idofn,ipoin) + solve(1) % block_array(1) % bvnat(idofn,ipoin)
                      end if
                   end do
                   do idofn = 1,ndof2
                      if( solve(1) % block_array(2) % kfl_fixno(idofn,ipoin) <= 0 ) then
                         b2(idofn,ipoin) = b2(idofn,ipoin) + solve(1) % block_array(2) % bvnat(idofn,ipoin)
                      end if
                   end do
                end do
             else
                do ipoin = 1,npoin
                   b1(1:ndof1,ipoin) = b1(1:ndof1,ipoin) + solve(1) % block_array(1) % bvnat(1:ndof1,ipoin)
                   b2(1:ndof2,ipoin) = b2(1:ndof2,ipoin) + solve(1) % block_array(2) % bvnat(1:ndof2,ipoin)
                end do
             end if
          end if

       case default
          !
          ! Do not know what to do
          !
          call runend('1BY1: DOES KNOW WHAT TO DO')

       end select

    end if

  end subroutine solver_reaction_force_2by2

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/10/2014
  !> @brief   Impose Dirichlet b.c.
  !> @details Impose Dirichlet b.c. on RHS and matrix
  !>
  !>          BEWARE - ALL THIS EXPLANATION HAS ERRORS - IT NEEDS TO BE CORRECTED ___ BETTER LOOK AT THE SUBROUTINE ITSELF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>          If SOLVE(1) % KFL_IFFIX = 1. For all IPOIN,
  !>
  !>            IF SOLVE(1) % KFL_FIXNO(:,IPOIN) > 0
  !>              - Residual-based solvers:
  !>                B1 (:,IPOIN) = 0
  !>              - Matrix-based solvers:
  !>                A11(:,:,IZDOM) = 0
  !>                A11(1,1,IZDOM) = 1
  !>                B1 (1,IPOIN)   = SOLVE(1) % BVESS(:,IPOIN)
  !>                U1 (1,IPOIN)   = SOLVE(1) % BVESS(:,IPOIN)
  !>            END IF
  !>
  !>          ELSE IF ( SOLVE(1) % KFL_IFFIX = 1 .and. .not. associated(SOLVE(1) % BVESS )
  !>
  !>            IF SOLVE(1) % KFL_FIXNO(:,IPOIN) > 0
  !>                A11(:,:,IZDOM) = 0
  !>                A11(1,1,IZDOM) = 1
  !>                B1 (1,IPOIN)   = SOLVE(1) % BVESS(:,IPOIN)   !hhh this does not make sense if it is not associated !!!!!!
  !>            END IF
  !>
  !>          ELSE If SOLVE(1) % KFL_IFFIX = 2
  !>
  !>            IF SOLVE(1) % KFL_FIXNO(:,IPOIN) > 0
  !>                A11(:,:,IZDOM) = 0
  !>                A11(1,1,IZDOM) = 1
  !>            END IF
  !>
  !>          END IF
  !>
  !>          +-         -+ +-  -+   +-  -+
  !>          |  Aii Aib  | | xi |   | bi |
  !>          |           | |    | = |    |
  !>          |  Abi Abb  | | xb |   | bb |
  !>          +-         -+ +-  -+   +-  -+
  !>
  !>          While the Dirichlet boundary condition is xb=xb', with
  !>
  !>          KFL_IFFIX = 1 : BVESS
  !>          KFL_IFFIX = 2 : xb = 0
  !>
  !>          To prescribed Dirichlet boudnary conditions, we have two options: 
  !>
  !>          KFL_DIRICHLET = 1 : Eliminate and enforce Dirichlet condition
  !>
  !>          +-         -+ +-  -+   +-            -+
  !>          |  Aii  0   | | xi |   | bi - Aib*xb' |
  !>          |           | |    | = |              |
  !>          |   0   D   | | xb |   | D*xb'        |
  !>          +-         -+ +-  -+   +-            -+
  !>
  !>          KFL_DIRICHLET = 2 : Eliminate Dirichlet condition
  !>
  !>          +-         -+ +-  -+   +-            -+
  !>          |  Aii  0   | | xi |   | bi - Aib*xb' |
  !>          |           | |    | = |              |
  !>          |  Aib  0   | | xb |   | bb - Abb*xb' |
  !>          +-         -+ +-  -+   +-            -+
  !>
  !>          The second option has the advantage that the residual of
  !>          the xb equation is the reaction force itself. The main drawback
  !>          is that a direct solver cannot be used without modification
  !>          as the matrix is singular.
  !>
  !>          This option should be used when
  !>          considering implicit coupling in order to comptue correctly the
  !>          residual on the interface when an interface node is a Dirichlet
  !>          boundary condition.
  !>
  !----------------------------------------------------------------------

  subroutine solver_preprocess_dirichlet_1by1(wtask,solve,ndof1,A11,b1,u1)
    character(*),      intent(in)             :: wtask
    type(soltyp),      intent(in),    pointer :: solve(:)
    integer(ip),       intent(in)             :: ndof1
    real(rp),          intent(inout)          :: A11(ndof1,ndof1,*)
    real(rp),          intent(inout)          :: b1(ndof1,*)
    real(rp),          intent(inout)          :: u1(ndof1,*)
    integer(ip)                               :: lpoin,ipoin,idofn,jdofn
    integer(ip)                               :: izdod,izdom,jpoin,kpoin
    integer(ip)                               :: jzdom,kfl_dirichlet
    real(rp)                                  :: amatrd(ndof1)
    integer(ip),                      pointer :: kfl_fixno(:,:)
    real(rp),                         pointer :: bvess(:,:)

    nullify(kfl_fixno)
    nullify(bvess)
    if( solve(1) % kfl_iffix == 0 ) return
    kfl_fixno     => solve(1) % kfl_fixno
    bvess         => solve(1) % bvess
    kfl_dirichlet =  solve(1) % kfl_dirichlet

    select case ( trim(wtask) )

    case ( 'ONLY RHS' )
       !
       ! Richardson residual-based solvers: put RHS=0 on Dirichlet nodes
       !
       do ipoin = 1,npoin
          do idofn = 1,ndof1
             if( kfl_fixno(idofn,ipoin) > 0 ) then
                b1(idofn,ipoin) = 0.0_rp
             end if
          end do
       end do

    case ( 'MATRIX AND RHS' )
       !
       ! Matrix and RHS
       !
       if( solve(1) % kfl_iffix == 1 .and. associated(bvess) ) then
          !
          ! Dirichlet value is given in BVESS
          !
          do ipoin = 1,npoin

             do idofn = 1,ndof1

                if( kfl_fixno(idofn,ipoin) > 0 ) then

                   if( kfl_dirichlet == 1 ) then
                      !
                      ! Eliminate dof of IPOIN from other equations (JPOIN)
                      ! Keep rows unchanged in order to compute the reaction force
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         jpoin = c_dom(izdom)
                         if( ipoin /= jpoin ) then
                            do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                               kpoin = c_dom(jzdom)
                               if( kpoin == ipoin ) then
                                  do jdofn = 1,ndof1
                                     b1(jdofn,jpoin)        = b1(jdofn,jpoin) - A11(idofn,jdofn,jzdom) * bvess(idofn,ipoin)
                                     A11(idofn,jdofn,jzdom) = 0.0_rp
                                  end do
                               end if
                            end do
                         end if
                      end do
                      !
                      ! IZDOD: Diagonal
                      !
                      izdod = r_dom(ipoin) - 1
                      jpoin = 0
                      do while( jpoin /= ipoin )
                         izdod = izdod + 1
                         jpoin = c_dom(izdod)
                      end do
                      amatrd(idofn) = A11(idofn,idofn,izdod)
                      if( abs(amatrd(idofn)) < zeror ) amatrd(idofn) = 1.0_rp
                      !
                      ! Set line to zero
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         do jdofn = 1,ndof1
                            A11(jdofn,idofn,izdom) = 0.0_rp
                         end do
                      end do
                      !
                      ! Presrcibe value
                      !
                      A11(idofn,idofn,izdod) = amatrd(idofn)
                      b1(idofn,ipoin)        = bvess(idofn,ipoin) * amatrd(idofn)
                      u1(idofn,ipoin)        = bvess(idofn,ipoin)

                   else if( kfl_dirichlet == 2 ) then
                      !
                      ! Eliminate dof of IPOIN from other equations (JPOIN)
                      ! Keep rows unchanged in order to compute the reaction force
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         jpoin = c_dom(izdom)
                         do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                            kpoin = c_dom(jzdom)
                            if( kpoin == ipoin ) then
                               do jdofn = 1,ndof1
                                  b1(jdofn,jpoin)        = b1(jdofn,jpoin) - A11(idofn,jdofn,jzdom) * bvess(idofn,ipoin)
                                  A11(idofn,jdofn,jzdom) = 0.0_rp
                               end do
                            end if
                         end do
                      end do
                      !
                      ! Prescribe value
                      !                      
                      u1(idofn,ipoin) = bvess(idofn,ipoin)

                   end if

                end if

             end do

          end do

       else if( solve(1) % kfl_iffix == 2 .or. ( solve(1) % kfl_iffix == 1 .and. .not. associated(bvess) ) ) then
          !
          ! Dirichlet value is not given: Impose zero
          !
          do ipoin = 1,npoin

             do idofn = 1,ndof1

                if( kfl_fixno(idofn,ipoin) > 0 ) then
                   
                   if( kfl_dirichlet == 1 ) then
                      !
                      ! Eliminate dof of IPOIN from other equations (JPOIN)
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         jpoin = c_dom(izdom)
                         if( ipoin /= jpoin ) then
                            do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                               kpoin = c_dom(jzdom)
                               if( kpoin == ipoin ) then
                                  do jdofn = 1,ndof1
                                     A11(idofn,jdofn,jzdom) = 0.0_rp
                                  end do
                               end if
                            end do
                            
                         end if
                      end do
                      !
                      ! IZDOD: Diagonal
                      !
                      izdod = r_dom(ipoin) - 1
                      jpoin = 0
                      do while( jpoin /= ipoin )
                         izdod = izdod + 1
                         jpoin = c_dom(izdod)
                      end do
                      amatrd(idofn) = A11(idofn,idofn,izdod)
                      if( abs(amatrd(idofn)) < zeror ) amatrd(idofn) = 1.0_rp
                      !
                      ! Set line to zero
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         do jdofn = 1,ndof1
                            A11(jdofn,idofn,izdom) = 0.0_rp
                         end do
                      end do
                      !
                      ! Prescribe value
                      !
                      A11(idofn,idofn,izdod) = amatrd(idofn)
                      b1(idofn,ipoin)        = 0.0_rp
                      u1(idofn,ipoin)        = 0.0_rp
                      
                   else if( kfl_dirichlet == 2 ) then
                      !
                      ! Eliminate dof of IPOIN from other equations (JPOIN)
                      ! Keep rows unchanged in order to compute the reaction force
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         jpoin = c_dom(izdom)
                         do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                            kpoin = c_dom(jzdom)
                            if( kpoin == ipoin ) then
                               do jdofn = 1,ndof1
                                  A11(idofn,jdofn,jzdom) = 0.0_rp
                               end do
                            end if
                         end do
                      end do
                      !
                      ! Prescribe value
                      !                      
                      u1(idofn,ipoin) = 0.0_rp
                      
                   end if
                   
                end if

             end do

          end do

       else if( solve(1) % kfl_iffix == 3 ) then
          !
          ! Dirichlet value is given in unknown u1
          !
          do ipoin = 1,npoin

             do idofn = 1,ndof1

                if( kfl_fixno(idofn,ipoin) > 0 ) then

                   if( kfl_dirichlet == 1 ) then                   
                      !
                      ! Eliminate dof of IPOIN from other equations (JPOIN)
                      ! Keep rows unchanged in order to compute the reaction force
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         jpoin = c_dom(izdom)
                         if( ipoin /= jpoin ) then
                            do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                               kpoin = c_dom(jzdom)
                               if( kpoin == ipoin ) then
                                  do jdofn = 1,ndof1
                                     b1(jdofn,jpoin)        = b1(jdofn,jpoin) - A11(idofn,jdofn,jzdom) * u1(idofn,ipoin)
                                     A11(idofn,jdofn,jzdom) = 0.0_rp
                                  end do
                               end if
                            end do
                         end if
                      end do
                      !
                      ! IZDOD: Diagonal
                      !
                      izdod = r_dom(ipoin) - 1
                      jpoin = 0
                      do while( jpoin /= ipoin )
                         izdod = izdod + 1
                         jpoin = c_dom(izdod)
                      end do
                      amatrd(idofn) = A11(idofn,idofn,izdod)
                      if( abs(amatrd(idofn)) < zeror ) amatrd(idofn) = 1.0_rp
                      !
                      ! Set line to zero
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         do jdofn = 1,ndof1
                            A11(jdofn,idofn,izdom) = 0.0_rp
                         end do
                      end do
                      !
                      ! Presrcibe value
                      !
                      A11(idofn,idofn,izdod) = amatrd(idofn)
                      b1(idofn,ipoin)        = u1(idofn,ipoin) * amatrd(idofn)

                   else if( kfl_dirichlet == 2 ) then
                      !
                      ! Eliminate dof of IPOIN from other equations (JPOIN)
                      ! Keep rows unchanged in order to compute the reaction force
                      !
                      do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                         jpoin = c_dom(izdom)
                         do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                            kpoin = c_dom(jzdom)
                            if( kpoin == ipoin ) then
                               do jdofn = 1,ndof1
                                  b1(jdofn,jpoin)        = b1(jdofn,jpoin) - A11(idofn,jdofn,jzdom) * u1(idofn,ipoin)
                                  A11(idofn,jdofn,jzdom) = 0.0_rp
                               end do
                            end if
                         end do
                      end do

                   end if
                   
                end if

             end do

          end do

       end if

    end select

  end subroutine solver_preprocess_dirichlet_1by1

!!$  !----------------------------------------------------------------------
!!$  !
!!$  !> @author  Guillaume Houzeaux
!!$  !> @date    13/04/2014
!!$  !> @brief   Impose Dirichlet boundary conditions
!!$  !> @details Impose Dirichlet boundary conditions on a 2x2 block
!!$  !>          matrix. Rotate if required.
!!$  !>
!!$  !----------------------------------------------------------------------
!!$
!!$  subroutine solver_preprocess_dirichlet_2by2(solve,ndof1,ndof2,A11,A12,A21,A22,P,b1,b2)
!!$    implicit none
!!$    type(soltyp), intent(in),  pointer  :: solve(:)
!!$    integer(ip),  intent(in)            :: ndof1
!!$    integer(ip),  intent(in)            :: ndof2
!!$    real(rp),     intent(out)           :: A11(ndof1,ndof1,nzdom)
!!$    real(rp),     intent(out)           :: A12(ndof1,nzdom)
!!$    real(rp),     intent(out)           :: A21(ndof1,nzdom)
!!$    real(rp),     intent(out)           :: A22(nzdom)
!!$    real(rp),     intent(out), optional :: P(nzdom)
!!$    real(rp),     intent(out)           :: b1(ndime,npoin)
!!$    real(rp),     intent(out)           :: b2(npoin)
!!$    real(rp)                            :: A11d(3),Qd,Appd
!!$    integer(ip)                         :: ipoin,jzdom,idime,jdime,izdom,jpoin
!!$    integer(ip)                         :: izdod,ibopo,jbopo,kpoin,iroty,kdime
!!$    integer(ip)                         :: kfl_addma,idofn
!!$    integer(4)                          :: istat
!!$    real(rp)                            :: worma(ndime,ndime)
!!$    real(rp),                  pointer  :: rotma(:,:)
!!$
!!$    if( ndof2 /= 1 ) then
!!$       call runend('PREPROCESS_DIRICHLET_2BY2: NOT CODED')
!!$    end if
!!$    if( ndof1 /= ndime .and. solve(1) % kfl_rotation == 1 ) then
!!$       call runend('PREPROCESS_DIRICHLET_2BY2: IMPOSSIBLE ROATION')
!!$    end if
!!$
!!$    if( solve(1) % kfl_iffix /= 0 ) then
!!$
!!$       !----------------------------------------------------------------------
!!$       !
!!$       ! Rotate matrix
!!$       !
!!$       !----------------------------------------------------------------------
!!$
!!$       if( solve(1) % kfl_rotation == 1 ) then
!!$
!!$          do ipoin = 1,npoin
!!$             ibopo = lpoty(ipoin)
!!$             if( ibopo > 0 ) then
!!$                iroty =  solve(1) % rotation_type(ibopo)
!!$                if( iroty /= 0 ) then
!!$                   rotma => solve(1) % rotation_matrix(:,:,iroty)
!!$                   !if( iroty == -1 ) then                                    ! Tangent system
!!$                   !   rotma => exnor(:,:,ibopo)
!!$                   !else if( iroty >= 1 ) then                                ! Given system
!!$                   !   rotma =>  skcos(:,:,iroty)
!!$                   !else if( iroty == -2 ) then                               ! Given system
!!$                   !   rotma => skcos_nsi(:,:,ibopo)
!!$                   !else if( iroty == -3 ) then                               ! Geometrical normal
!!$                   !   rotma => skcos(:,:,ibopo)
!!$                   !end if
!!$                   !
!!$                   ! Modifies column number IPOIN of AMATR ( A_j,imodi <-- A_j,imodi R )
!!$                   !
!!$                   do jpoin = 1,npoin
!!$                      do izdom = r_dom(jpoin),r_dom(jpoin+1)-1
!!$                         kpoin = c_dom(izdom)
!!$                         if( kpoin == ipoin ) then
!!$                            do idime = 1,ndime
!!$                               do jdime = 1,ndime
!!$                                  worma(idime,jdime) = 0.0_rp
!!$                                  do kdime = 1,ndime
!!$                                     worma(idime,jdime) = worma(idime,jdime) &
!!$                                          + A11(kdime,idime,izdom) * rotma(kdime,jdime)
!!$                                  end do
!!$                               end do
!!$                            end do
!!$                            do idime = 1,ndime
!!$                               do jdime = 1,ndime
!!$                                  A11(jdime,idime,izdom) = worma(idime,jdime)
!!$                               end do
!!$                            end do
!!$
!!$                            do jdime = 1,ndime
!!$                               worma(1,jdime) = 0.0_rp
!!$                               do kdime = 1,ndime
!!$                                  worma(1,jdime) = worma(1,jdime)&
!!$                                       + A21(kdime,izdom) * rotma(kdime,jdime)
!!$                               end do
!!$                            end do
!!$                            do jdime = 1,ndime
!!$                               A21(jdime,izdom) = worma(1,jdime)
!!$                            end do
!!$
!!$                         end if
!!$                      end do
!!$                   end do
!!$
!!$                   do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
!!$                      !
!!$                      ! Modifies row number IPOIN of AMATR ( A_imodi,j <-- R^t A_imodi,j )
!!$                      !
!!$                      jpoin = c_dom(izdom)
!!$                      do idime = 1,ndime
!!$                         do jdime = 1,ndime
!!$                            worma(idime,jdime) = 0.0_rp
!!$                            do kdime = 1,ndime
!!$                               worma(idime,jdime) = worma(idime,jdime) &
!!$                                    + A11(jdime,kdime,izdom) * rotma(kdime,idime)
!!$                            end do
!!$                         end do
!!$                      end do
!!$                      do idime = 1,ndime
!!$                         do jdime = 1,ndime
!!$                            A11(jdime,idime,izdom) = worma(idime,jdime)
!!$                         end do
!!$                      end do
!!$                      !
!!$                      ! Modify the part corresponding to a scalar unknown
!!$                      !
!!$                      do idime = 1,ndime
!!$                         worma(idime,1) = 0.0_rp
!!$                         do kdime = 1,ndime
!!$                            worma(idime,1) = worma(idime,1) &
!!$                                 + rotma(kdime,idime) * A12(kdime,izdom)
!!$                         end do
!!$                      end do
!!$                      do idime = 1,ndime
!!$                         A12(idime,izdom) = worma(idime,1)
!!$                      end do
!!$
!!$                   end do
!!$                   !
!!$                   ! Rotate RHS: bu
!!$                   !
!!$                   do idime = 1,ndime
!!$                      worma(idime,1) = 0.0_rp
!!$                      do kdime = 1,ndime
!!$                         worma(idime,1) = worma(idime,1) &
!!$                              + rotma(kdime,idime) * b1(kdime,ipoin)
!!$                      end do
!!$                   end do
!!$                   do idime = 1,ndime
!!$                      b1(idime,ipoin) = worma(idime,1)
!!$                   end do
!!$
!!$                end if
!!$             end if
!!$          end do
!!$       end if
!!$
!!$       !----------------------------------------------------------------------
!!$       !
!!$       ! Impose first degree of freedom
!!$       !
!!$       !----------------------------------------------------------------------
!!$
!!$       do ipoin = 1,npoin
!!$
!!$          do idof1 = 1,ndof1
!!$
!!$             if( solve_sol(1) % kfl_fixno(idof1,ipoin) > 0 ) then
!!$                !
!!$                ! Eliminate dof of IPOIN from other equations (JPOIN)
!!$                !
!!$                do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
!!$                   jpoin = c_dom(izdom)
!!$                   if( ipoin /= jpoin ) then
!!$
!!$                      do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
!!$                         kpoin = c_dom(jzdom)
!!$                         if( kpoin == ipoin ) then
!!$                            do jdof1 = 1,ndof1
!!$                               b1(jdof1,jpoin)        = b1(jdof1,jpoin) - A11(idof1,jdof1,jzdom) * solve_sol(1) % bvess(idof1,ipoin)
!!$                               A11(idof1,jdof1,jzdom) = 0.0_rp
!!$                            end do
!!$                            b2(jpoin)        = b2(jpoin) - A21(idof1,jzdom) * solve_sol(1) % bvess(idof1,ipoin)
!!$                            A21(idof1,jzdom) = 0.0_rp
!!$                         end if
!!$                      end do
!!$
!!$                   end if
!!$                end do
!!$                !
!!$                ! IZDOD: Diagonal
!!$                !
!!$                izdod = r_dom(ipoin) - 1
!!$                jpoin = 0
!!$                do while( jpoin /= ipoin )
!!$                   izdod = izdod + 1
!!$                   jpoin = c_dom(izdod)
!!$                end do
!!$                A11d(idof1) = A11(idof1,idof1,izdod)
!!$                if( abs(A11d(idof1)) < zeror ) A11d(idof1) = 1.0_rp
!!$                !
!!$                ! Set line to zero
!!$                !
!!$                do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
!!$                   do jdof1 = 1,ndof1
!!$                      A11(jdof1,idof1,izdom) = 0.0_rp
!!$                   end do
!!$                   A12(idof1,izdom) = 0.0_rp
!!$                end do
!!$                !
!!$                ! Prescribe value
!!$                !
!!$                idofn                  = (ipoin-1)*ndof1 + idof1
!!$                A11(idof1,idof1,izdod) = A11d(idof1)
!!$                b1(idof1,ipoin)        = solve_sol(1) % bvess(idof1,ipoin) * A11d(idof1)
!!$
!!$             end if
!!$
!!$          end do
!!$
!!$       end do
!!$
!!$       !----------------------------------------------------------------------
!!$       !
!!$       ! Impose Schur complement preconditioner for second system
!!$       !
!!$       !----------------------------------------------------------------------
!!$
!!$       if( present(P) ) then
!!$
!!$          if( solve(2) % kfl_symme == 1 ) then
!!$
!!$             do ipoin = 1,npoin
!!$                do izdom = r_sym(ipoin),r_sym(ipoin+1) - 2
!!$                   jpoin = c_sym(izdom)
!!$                   if( solve_sol(2) % kfl_fixno(1,jpoin) > 0 ) then
!!$                      P(izdom) = 0.0_rp
!!$                   end if
!!$                end do
!!$             end do
!!$             do ipoin = 1,npoin
!!$                if( solve_sol(2) % kfl_fixno(1,ipoin) > 0 ) then
!!$                   !
!!$                   ! IZDOD: Diagonal
!!$                   !
!!$                   izdod = r_sym(ipoin+1) - 1
!!$                   Qd = P(izdod)
!!$                   if( abs(Qd) < zeror ) Qd = 1.0_rp
!!$                   !
!!$                   ! Set line to zero
!!$                   !
!!$                   do izdom = r_sym(ipoin),r_sym(ipoin+1) - 1
!!$                      P(izdom) = 0.0_rp
!!$                   end do
!!$                   !
!!$                   ! Presrcibe value
!!$                   !
!!$                   P(izdod) = Qd
!!$                end if
!!$             end do
!!$
!!$          else
!!$
!!$             do ipoin = 1,npoin
!!$                do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
!!$                   jpoin = c_dom(izdom)
!!$                   if( ipoin /= jpoin ) then
!!$                      if( solve_sol(2) % kfl_fixno(1,jpoin) > 0 ) then
!!$                         P(izdom) = 0.0_rp
!!$                      end if
!!$                   end if
!!$                end do
!!$             end do
!!$
!!$             do ipoin = 1,npoin
!!$                if( solve_sol(2) % kfl_fixno(1,ipoin) > 0 ) then
!!$                   !
!!$                   ! IZDOD: Diagonal
!!$                   !
!!$                   izdod = r_dom(ipoin) - 1
!!$                   jpoin = 0
!!$                   do while( jpoin /= ipoin )
!!$                      izdod = izdod + 1
!!$                      jpoin = c_dom(izdod)
!!$                   end do
!!$                   Qd = P(izdod)
!!$                   if( abs(Qd) < zeror ) Qd = 1.0_rp
!!$                   !
!!$                   ! Set line to zero
!!$                   !
!!$                   do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
!!$                      P(izdom) = 0.0_rp
!!$                   end do
!!$                   !
!!$                   ! Presrcibe value
!!$                   !
!!$                   P(izdod) = Qd
!!$                end if
!!$             end do
!!$
!!$          end if
!!$
!!$       end if
!!$
!!$    end if
!!$
!!$    !----------------------------------------------------------------------
!!$    !
!!$    ! Impose pressure
!!$    !
!!$    !----------------------------------------------------------------------
!!$
!!$    if( kfl_matdi_nsi == 1 .and. INOTMASTER ) then
!!$
!!$       if ( kfl_confi_nsi == 1 .and. nodpr_nsi > 0 ) then
!!$
!!$          ipoin = nodpr_nsi
!!$
!!$          if( solve(2) % kfl_symme == 1 ) then
!!$             call runend('NOT CODED: CHECK IT')
!!$             !
!!$             ! Eliminate pressure at IPOIN on all lines
!!$             !
!!$             do jpoin = 1,npoin
!!$                if( ipoin /= jpoin ) then
!!$                   izdom = r_sym(jpoin)
!!$                   do while( izdom < r_sym(jpoin+1) )
!!$                      kpoin = c_sym(izdom)
!!$                      if( kpoin == ipoin ) then
!!$                         b2(jpoin)  = b2(jpoin) - A22(izdom) * valpr_nsi
!!$                         A22(izdom) = 0.0_rp
!!$                         P(izdom)   = 0.0_rp
!!$                         izdom      = r_sym(jpoin+1)
!!$                      end if
!!$                      izdom = izdom + 1
!!$                   end do
!!$                end if
!!$             end do
!!$             !
!!$             ! Diagonal
!!$             !
!!$             izdod = r_sym(ipoin+1) - 1
!!$             Qd   = P(izdod)
!!$             Appd = A22(izdod)
!!$             if( abs(Qd)   < zeror )   Qd = 1.0_rp
!!$             if( abs(Appd) < zeror ) Appd = 1.0_rp
!!$             !
!!$             ! Set line to zero
!!$             !
!!$             do izdom = r_sym(ipoin),r_sym(ipoin+1) - 1
!!$                P(izdom)   = 0.0_rp
!!$                A22(izdom) = 0.0_rp
!!$             end do
!!$             do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
!!$                do idof1 = 1,ndof1
!!$                   A21(idof1,izdom) = 0.0_rp
!!$                end do
!!$             end do
!!$             !
!!$             ! Presrcibe value
!!$             !
!!$             P(izdod)   = Qd
!!$             A22(izdod) = Appd
!!$             b2(ipoin)  = Appd * valpr_nsi
!!$
!!$          else
!!$             !
!!$             ! Eliminate pressure at IPOIN on all lines
!!$             !
!!$             do jpoin = 1,npoin
!!$                if( ipoin /= jpoin ) then
!!$                   izdom = r_dom(jpoin)
!!$                   do while( izdom < r_dom(jpoin+1) )
!!$                      kpoin = c_dom(izdom)
!!$                      if( kpoin == ipoin ) then
!!$                         b2(jpoin)  = b2(jpoin) - A22(izdom) * valpr_nsi
!!$                         A22(izdom) = 0.0_rp
!!$                         P(izdom)   = 0.0_rp
!!$                         izdom      = r_dom(jpoin+1)
!!$                      end if
!!$                      izdom = izdom + 1
!!$                   end do
!!$                end if
!!$             end do
!!$             !
!!$             ! Diagonal
!!$             !
!!$             izdod = r_dom(ipoin) - 1
!!$             jpoin = 0
!!$             do while( jpoin /= ipoin )
!!$                izdod = izdod + 1
!!$                jpoin = c_dom(izdod)
!!$             end do
!!$             Appd = A22(izdod)
!!$             Qd   = P(izdod)
!!$             if( abs(Appd) < zeror ) Appd = 1.0_rp
!!$             !
!!$             ! Set line to zero
!!$             !
!!$             do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
!!$                A22(izdom) = 0.0_rp
!!$                P(izdom)   = 0.0_rp
!!$                do idof1 = 1,ndof1
!!$                   A21(idof1,izdom) = 0.0_rp
!!$                end do
!!$             end do
!!$             !
!!$             ! Presrcibe value
!!$             !
!!$             P(izdod)   = Qd
!!$             A22(izdod) = Appd

!!$             b2(ipoin)  = Appd * valpr_nsi
!!$
!!$          end if
!!$          !
!!$          ! Eliminate prescribed pressure from momentum equation
!!$          !
!!$          do jpoin = 1,npoin
!!$             izdom = r_dom(jpoin)
!!$             do while( izdom < r_dom(jpoin+1) )
!!$                kpoin = c_dom(izdom)
!!$                if( kpoin == ipoin ) then
!!$                   do idof1 = 1,ndof1
!!$                      b1(idof1,jpoin)  = b1(idof1,jpoin) - A12(idof1,izdom) * valpr_nsi
!!$                      A12(idof1,izdom) = 0.0_rp
!!$                   end do
!!$                   izdom = r_dom(jpoin+1)
!!$                end if
!!$                izdom = izdom + 1
!!$             end do
!!$          end do
!!$
!!$       end if
!!$
!!$    end if
!!$
!!$  end subroutine solver_preprocess_dirichlet_2by2

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    23/05/2017
  !> @brief   Allocate memory
  !> @details Allocate memory for an algebraic system matrix A and
  !>          possibly right-hand side b
  !>
  !----------------------------------------------------------------------

  subroutine solver_allocate_system_memory(solve,an,bb)

    type(soltyp), intent(in)                       :: solve
    real(rp),     intent(inout), pointer           :: an(:)
    real(rp),     intent(inout), pointer, optional :: bb(:)
    integer(ip)                                    :: nnz,nn,ndofn

    nn    = solve % nunkn
    nnz   = solve % nnz
    ndofn = solve % ndofn

    call memory_alloca(memit,'AN','solver_allocate_system_memory',an,nnz*ndofn*ndofn)
    if( present(bb) ) then
       call memory_alloca(memit,'AN','solver_allocate_system_memory',bb,nn)
    end if

  end subroutine solver_allocate_system_memory

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    23/05/2017
  !> @brief   Diagonal of a matrix
  !> @details Compute the diagonal of a matrix
  !>
  !----------------------------------------------------------------------

  subroutine solver_diagonal(solve,an,diagonal)

    type(soltyp), intent(in)  :: solve(:)
    real(rp),     intent(in)  :: an(*)
    real(rp),     intent(out) :: diagonal(*)
    integer(ip)               :: nn,ndof

    nn   = solve(1) % nequa
    ndof = solve(1) % ndofn
       
    if( INOTMASTER .and. nn > 0 ) then
       
       if(      solve(1) % kfl_format == SOL_CSR_FORMAT ) then
          !
          ! CSR format
          !
          if( solve(1) % kfl_full_rows == 0 ) then
             call matrix_diagonal_CSR(nn,ndof,solve(1) % kfl_symme,solve(1) % ia,solve(1) % ja,an,diagonal)
          else
             call matrix_diagonal_CSR(nn,ndof,solve(1) % kfl_symme,solve(1) % ia_full,solve(1) % ja_full,an,diagonal)
          end if

       else if( solve(1) % kfl_format == SOL_COO_FORMAT ) then
          !
          ! COO format
          !
          call matrix_diagonal_COO(nn,ndof,solve(1) % rows,solve(1) % cols,an,diagonal)

       else if( solve(1) % kfl_format == SOL_ELL_FORMAT ) then
          !
          ! ELL format
          !
          call matrix_diagonal_ELL(nn,ndof,solve(1) % cols_ell,an,diagonal)

       end if
       !
       ! Periodicity and Parallelization
       !
       if( solve(1) % kfl_full_rows == 0 ) call pararr('SOL',NPOIN_TYPE,ndof*nn,diagonal)

    end if

  end subroutine solver_diagonal

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    20/11/2017
  !> @brief   Assemble element matrix
  !> @details Assemble an element matrix ELMAT into the global matrix
  !>          AN for different format: CSR, COO, ELL
  !>
  !----------------------------------------------------------------------

  subroutine solver_assemble_element_matrix(&
       solve,ndof,pnode,pevat,ielem,lnods,elmat,an)

    use def_kermod, only : kfl_element_to_csr
    use def_domain, only : lezdo

    type(soltyp), intent(in)             :: solve(:)
    integer(ip),  intent(in)             :: ndof
    integer(ip),  intent(in)             :: pnode
    integer(ip),  intent(in)             :: pevat
    integer(ip),  intent(in)             :: ielem
    integer(ip),  intent(in)             :: lnods(pnode)
    real(rp),     intent(in)             :: elmat(pevat,pevat)
    real(rp),     intent(inout)          :: an(ndof,ndof,*)

    if(      solve(1) % kfl_format == SOL_CSR_FORMAT ) then
       !
       ! CSR formats
       !
       call  matrix_assemble_element_matrix_to_CSR(&
            kfl_element_to_csr,ndof,pnode,pevat,&
            ielem,lnods,elmat,solve(1) % ia,solve(1) % ja,an,lezdo)

    else if( solve(1) % kfl_format == SOL_COO_FORMAT ) then
       !
       ! COO format
       !
       call  matrix_assemble_element_matrix_to_COO(&
            kfl_element_to_csr,ndof,pnode,pevat,&
            ielem,lnods,elmat,solve(1) % rows,solve(1) % cols,an,lezdo)

    else if( solve(1) % kfl_format == SOL_ELL_FORMAT ) then
       !
       ! ELL format
       !
       call matrix_assemble_element_matrix_to_ELL(&
            ndof,pnode,pevat,lnods,elmat,solve(1) % cols_ell,an)

    end if

  end subroutine solver_assemble_element_matrix

  !------------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    22/11/2017
  !> @brief   Parallel SpMV
  !> @details Parallel SpMV
  !>          Multiply a non symmetric matrix by a vector YY = A XX
  !>          where A is in CSR, COO or ELL format. This is the parallel
  !>          version, including blocking and non-blocking send-receive
  !>          MPI functions.
  !>
  !>          INPUT
  !>             NBNODES .... Number of equations
  !>             NBVAR ...... Number of variables
  !>             AN ......... Matrix
  !>             JA ......... List of elements
  !>             IA ......... Pointer to list of elements
  !>             XX ......... Vector
  !>          OUTPUT
  !>             YY ......... result vector
  !>
  !>          Full row:
  !>          ---------
  !>
  !>          +-     -+   +-                    -+ +-     -+
  !>          | y_own |   | A_own_own  A_own_oth | | x_own |
  !>          |       | = |                      | |       |
  !>          | y_oth |   | A_oth_own  A_oth_oth | | x_oth |
  !>          +-     -+   +-                    -+ +-     -+
  !>
  !>          1. Receive x_oth from neighbors
  !>          2. y_own^* = A_own_own x_own
  !>          3. y_oth^* = A_oth_own x_own
  !>          4. Waitall
  !>          5. y_own   = y_own^* + A_own_oth x_oth
  !>          6. y_oth   = y_oth^* + A_oth_oth x_oth
  !>
  !> @}
  !------------------------------------------------------------------------

  subroutine solver_parallel_SpMV(&
       solve,an,xx,yy,&
       MPI,OPENMP,INITIALIZATION,TIMING,MY_TIMING,ASYNCHRONISM,OPENMP_INTERFACE)

    type(soltyp), intent(inout)         :: solve            !< Solver type
    real(rp),     intent(in)            :: an(*)            !< Matrix
    real(rp),     intent(inout)         :: xx(*)            !< Multiplicand
    real(rp),     intent(out)           :: yy(*)            !< Result vector
    logical(lg),  intent(in),  optional :: MPI              !< If MPI should be used or not
    logical(lg),  intent(in),  optional :: OPENMP           !< If OpenMP should be used or not
    logical(lg),  intent(in),  optional :: INITIALIZATION   !< If result vector should be initialized
    logical(lg),  intent(in),  optional :: TIMING           !< If timing should be activated 
    real(rp),     intent(out), optional :: MY_TIMING(2)     !< Register timing in this array
    logical(lg),  intent(in),  optional :: ASYNCHRONISM     !< Blocking or non-blocking communications
    logical(lg),  intent(in),  optional :: OPENMP_INTERFACE !< If OpenMP should be used or not on interface
    integer(ip)                         :: CHUNK          
    integer(ip)                         :: CHUNK_INTERFACE          
    character(7)                        :: SCHEDULE       

    integer(ip)                         :: nn,ndof

    real(rp)                            :: cpu_computation
    real(rp)                            :: cpu_communication
    real(rp)                            :: time1,time2,time3
    real(rp)                            :: time4,time5

    logical(lg)                         :: do_timing
    logical(lg)                         :: do_openmp
    logical(lg)                         :: do_openmp_interface
    logical(lg)                         :: do_init
    logical(lg)                         :: do_mpi
    logical(lg)                         :: do_asynchronism

    nn   = solve % nequa
    ndof = solve % ndofn
       
    if( INOTMASTER .and. nn > 0 ) then

       !-------------------------------------------------------------------
       !
       ! Options
       !
       !-------------------------------------------------------------------

       do_mpi              = .true.
       do_openmp           = .false.
       do_init             = .true.
       do_timing           = .true.
       if( kfl_async == 1 ) then
          do_asynchronism = .true.
       else
          do_asynchronism = .false.
       end if

       if( present(MPI)              ) do_mpi          = MPI       
       if( present(OPENMP)           ) do_openmp       = OPENMP
       if( present(TIMING)           ) do_timing       = TIMING
       if( present(MY_TIMING)        ) do_timing       = .true.      
       if( present(INITIALIZATION)   ) do_init         = INITIALIZATION
       if( present(ASYNCHRONISM)     ) do_asynchronism = ASYNCHRONISM
       !
       ! OpenMP
       !
       if( solve % omp_schedule == SOL_OMP_OFF ) then
          do_openmp = .false.
       else
          do_openmp = .true.
       end if
       if( do_openmp ) then
          if(      solve % omp_schedule == SOL_OMP_STATIC ) then
             SCHEDULE = 'STATIC'
          else if( solve % omp_schedule == SOL_OMP_GUIDED ) then
             SCHEDULE = 'GUIDED'
          else if( solve % omp_schedule == SOL_OMP_DYNAMIC ) then
             SCHEDULE = 'DYNAMIC'
             CHUNK = solve % omp_chunk_size
          end if
       end if
       !
       ! OpenMP on interface
       !
       do_openmp_interface = do_openmp
       if( present(OPENMP_INTERFACE) ) do_openmp_interface = OPENMP_INTERFACE  
       if( do_openmp_interface ) then
          CHUNK_INTERFACE = 100          
       else
          CHUNK_INTERFACE = solve % omp_interface
       end if
       
       if( .not. do_mpi .or. ISEQUEN ) then

          !-------------------------------------------------------------------
          !
          ! Sequential or local SpMV
          !
          !-------------------------------------------------------------------

          call cputim(time1)

          if( solve % kfl_format == SOL_CSR_FORMAT ) then
             !
             ! CSR format
             !
             call matrix_CSR_SpMV(&
                  1_ip,nn,ndof,ndof,solve % ia,solve % ja,an,xx,yy,&
                  OPENMP=do_openmp,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)
                   
          else if( solve % kfl_format == SOL_COO_FORMAT ) then
             !
             ! COO format
             !
             call matrix_COO_SpMV(&
                  1_ip,nn,ndof,ndof,solve % rows,solve % cols,an,xx,yy,&
                  OPENMP=do_openmp,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)          

          else if( solve % kfl_format == SOL_ELL_FORMAT ) then
             !
             ! ELL format
             !
             call matrix_ELL_SpMV(&
                  1_ip,nn,ndof,ndof,solve % cols_ell,an,xx,yy,&
                  OPENMP=do_openmp,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)

          end if

          call cputim(time2)
          cpu_computation   = time2 - time1
          cpu_communication = 0.0_rp

       else

          !-------------------------------------------------------------------
          !
          ! Parallel
          !
          !-------------------------------------------------------------------
          !
          ! Exchange solution
          !
          if( solve % kfl_where == SOL_ELEMENTS ) then
             call PAR_GHOST_ELEMENT_EXCHANGE(ndof,xx,'SUM','IN MY CODE')
          end if

          if( solve % kfl_format == SOL_CSR_FORMAT ) then

             if( do_asynchronism ) then

                !-------------------------------------------------------------------
                !
                ! CSR format - Asynchronous
                !
                !-------------------------------------------------------------------

                if( solve % kfl_full_rows == 1 ) then
                   !
                   ! Full row
                   !
                   call cputim(time1)
                   !call PAR_INTERFACE_OWN_NODE_EXCHANGE(ndof,xx)

                   call PAR_INTERFACE_OWN_NODE_EXCHANGE(ndof,xx,'SEND RECEIVE')

                   call cputim(time2)
                   call matrix_CSR_SpMV(&
                        1_ip,nn,ndof,ndof,solve % ia_full,&
                        solve % ja_full,an,xx,yy,ia_opt=solve % ia_full_end,&
                        OPENMP=do_openmp,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)

                   call cputim(time3)
                   call PAR_INTERFACE_OWN_NODE_EXCHANGE(ndof,xx,'WAIT AND ASSEMBLE')

                   call cputim(time4)
                   call matrix_CSR_SpMV(&
                        1_ip,nn,ndof,ndof,solve % ia_full_ini,&
                        solve % ja_full,an,xx,yy,ia_opt=solve % ia_full,&
                        OPENMP=do_openmp_interface,INITIALIZATION=.false.,CHUNK=CHUNK,SCHEDULE=SCHEDULE)

                   call cputim(time5)
                   cpu_computation   = ( time3 + time5 ) - ( time2 + time4 )
                   cpu_communication = ( time2 + time4 ) - ( time1 + time3 )

                else
                   !
                   ! Partial row
                   !
                   call cputim(time1)
                   call matrix_CSR_SpMV(&
                        npoi1+1_ip,nn,ndof,ndof,solve % ia,solve % ja,an,xx,yy,&
                        OPENMP=do_openmp_interface,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)

                   call cputim(time2)
                   call pararr('SLA',NPOIN_TYPE,nn*ndof,yy)

                   call cputim(time3)
                   call matrix_CSR_SpMV(&
                        1_ip,npoi1,ndof,ndof,solve % ia,solve % ja,an,xx,yy,&
                        OPENMP=do_openmp,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)

                   call cputim(time4)
                   call pararr('SLA',NPOIN_TYPE,nn*ndof,yy)

                   call cputim(time5)

                   cpu_computation   = ( time2 + time4 ) - ( time1 + time3 )
                   cpu_communication = ( time5 + time3 ) - ( time4 + time2 )

                end if

             else

                !-------------------------------------------------------------------
                !
                ! CSR format - Synchronous
                !
                !-------------------------------------------------------------------

                if( solve % kfl_full_rows == 1 ) then
                   !
                   ! Full row
                   !
                   call cputim(time1)
                   call PAR_INTERFACE_OWN_NODE_EXCHANGE(ndof,xx)

                   call cputim(time2)
                   call matrix_CSR_SpMV(1_ip,nn,ndof,ndof,solve % ia_full,solve % ja_full,an,xx,yy,&
                        OPENMP=do_openmp,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)

                   call cputim(time3)
                   cpu_computation   = time3 - time2
                   cpu_communication = time2 - time1

                else
                   !
                   ! Partial row
                   !
                   call cputim(time1)
                   call matrix_CSR_SpMV(1_ip,nn,ndof,ndof,solve % ia,solve % ja,an,xx,yy,&
                        OPENMP=do_openmp,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)
                   !
                   ! Modify YY due do periodicity and Parall service
                   !
                   call cputim(time2)
                   if(      solve % kfl_where == SOL_NODES ) then
                      call pararr('SOL',NPOIN_TYPE,ndof*nn,yy)
                   else if( solve % kfl_where == SOL_EDGES ) then
                      call PAR_INTERFACE_EDGE_EXCHANGE(ndof,yy,'SUM','IN MY CODE')
                   end if
                   call cputim(time3)
                   cpu_computation   = time2 - time1
                   cpu_communication = time3 - time2

                end if

             end if

          else if( solve % kfl_format == SOL_COO_FORMAT ) then

             if( do_asynchronism ) then

                !-------------------------------------------------------------------
                !
                ! COO format: asynchronous
                !
                !-------------------------------------------------------------------

                call cputim(time1)
                call matrix_COO_SpMV(&
                     npoi1+1_ip,nn,ndof,ndof,solve % rows,solve % cols,an,xx,yy,&
                     OPENMP=do_openmp_interface,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)

                call cputim(time2)
                call pararr('SLA',NPOIN_TYPE,nn*ndof,yy)

                call cputim(time3)
                call matrix_COO_SpMV(&
                     1_ip,npoi1,ndof,ndof,solve % rows,solve % cols,an,xx,yy,&
                     OPENMP=do_openmp,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)

                call cputim(time4)
                call pararr('SLA',NPOIN_TYPE,nn*ndof,yy)

                call cputim(time5)
                cpu_computation   = ( time2 + time4 ) - ( time1 + time3 )
                cpu_communication = ( time5 + time3 ) - ( time4 + time2 )

             else

                !-------------------------------------------------------------------
                !
                ! COO format: synchronous
                !
                !-------------------------------------------------------------------

                call cputim(time1)
                call matrix_COO_SpMV(1_ip,nn,ndof,ndof,solve % rows,solve % cols,an,xx,yy,&
                     OPENMP=do_openmp,INITIALIZATION=do_init,CHUNK=CHUNK,SCHEDULE=SCHEDULE)
                !
                ! Modify YY due do periodicity and Parall service
                !
                call cputim(time2)
                if(      solve % kfl_where == SOL_NODES ) then
                   call pararr('SOL',NPOIN_TYPE,ndof*nn,yy)
                else if( solve % kfl_where == SOL_EDGES ) then
                   call PAR_INTERFACE_EDGE_EXCHANGE(ndof,yy,'SUM','IN MY CODE')
                end if
                call cputim(time3)
                cpu_computation   = time2 - time1
                cpu_communication = time3 - time2               

             end if

          else if( solve % kfl_format == SOL_ELL_FORMAT ) then

             if( do_asynchronism ) then

                !-------------------------------------------------------------------
                !
                ! ELL format: Asynchronous
                !
                !-------------------------------------------------------------------

                call cputim(time1)
                call matrix_ELL_SpMV(&
                     npoi1+1,nn,ndof,ndof,solve % cols_ell,an,xx,yy,&
                     OPENMP=do_openmp_interface,CHUNK=CHUNK)

                call cputim(time2)
                call pararr('SLA',NPOIN_TYPE,nn*ndof,yy)

                call cputim(time3)
                call matrix_ELL_SpMV(1_ip,npoi1,ndof,ndof,solve % cols_ell,an,xx,yy,&
                     OPENMP=do_openmp,CHUNK=CHUNK)

                call cputim(time4)
                call pararr('SLA',NPOIN_TYPE,nn*ndof,yy)

                call cputim(time5)
                cpu_computation   = ( time2 + time4 ) - ( time1 + time3 )
                cpu_communication = ( time5 + time3 ) - ( time4 + time2 )             

             else

                !-------------------------------------------------------------------
                !
                ! ELL format: Synchronous
                !
                !-------------------------------------------------------------------

                call cputim(time1)
                call matrix_ELL_SpMV(1_ip,nn,ndof,ndof,solve % cols_ell,an,xx,yy,&
                     OPENMP=do_openmp,CHUNK=CHUNK)

                call cputim(time2)
                call pararr('SOL',NPOIN_TYPE,ndof*nn,yy)

                call cputim(time3)
                cpu_computation   = time2 - time1
                cpu_communication = time3 - time2

             end if

          end if

       end if

       if( do_timing ) then
          if( present(MY_TIMING) ) then
             MY_TIMING(1)        = cpu_computation 
             MY_TIMING(2)        = cpu_communication
          else
             solve % cpu_spmv(1) = solve % cpu_spmv(1) + cpu_computation 
             solve % cpu_spmv(2) = solve % cpu_spmv(2) + cpu_communication
             solve % cpu_spmv(3) = solve % cpu_spmv(1) + solve % cpu_spmv(2)
             solve % num_spmv    = solve % num_spmv    + 1
          end if
       end if

    else if( .not. present(MY_TIMING) ) then

       if( present(MY_TIMING) ) then
          MY_TIMING(1)        = 0.0_rp
          MY_TIMING(2)        = 0.0_rp
       else
          solve % num_spmv    = solve % num_spmv    + 1
       end if
       
    end if

  end subroutine solver_parallel_SpMV

  !------------------------------------------------------------------------
  !> 
  !> @author  Guillaume Houzeaux
  !> @date    21/11/2017
  !> @brief   Conditioning
  !> @details Compute the conditioning of a matrix using the power method
  !>          also referred to as Von Mises iteration
  !>          \verbatim
  !>            y        = A x^k
  !>            x^k+1    = y/||y||
  !>            lambda^k = (x^k,y^k)/||x^k||
  !>          \endverbatim
  !>          If this algorith does not converge, return kappa=-1
  !>
  !>          -k*Lapl(u) = f en [0,1] ; h=1/(N+1); h=0.1
  !>
  !>          In finite element: 
  !>
  !>          lambda_i = 2.0*k/h*(1-cos(pi*i*h))
  !>          for k=1,N=9
  !>          lambda_min = 0.9788696741
  !>          lambda_max = 39.0211303259
  !>
  !>          In finite difference: 
  !>
  !>          lambda_i = 2.0*k/h**2*(1-cos(pi*i*h))
  !>
  !>          Example: 1D Laplacian
  !>          o---o---o---o---o---o---o---o---o---o---o
  !>          1   2   3   4   5   6   7   8   9  10  11
  !>          N=9
  !>          lambda_i   = 2[1-cos(pi*i/(N+1))]
  !>          lambda_min = 0.09788696741
  !>          lambda_max = 3.90211303259
  !>          A value of -888 means that the algorithm was not converged
  !> @}
  !------------------------------------------------------------------------

  subroutine solver_condition_number(solve,an,kappa,lambda_min,lambda_max)

    use def_master, only : INOTSLAVE
    
    type(soltyp), intent(inout)          :: solve      !< Solver structure
    real(rp),     intent(in)             :: an(*)      !< Matrix
    real(rp),     intent(out)            :: kappa      !> Condition number
    real(rp),     intent(out),  optional :: lambda_min !> Minimum eigenvalue
    real(rp),     intent(out),  optional :: lambda_max !> Maximum eigenvalue

    integer(ip)                          :: iiter
    integer(ip)                          :: ndof,nn
    integer(ip)                          :: nrows,ncols
    integer(ip)                          :: nsize,ii
    integer(ip),  parameter              :: maxit = 100
    real(rp)                             :: toler,eps
    real(rp),     pointer                :: xx(:)
    real(rp),     pointer                :: yy(:)
    real(rp)                             :: ynorm
    real(rp)                             :: xnorm
    real(rp)                             :: xnorm_new
    real(rp)                             :: xdoty
    real(rp)                             :: lambda,alpha
    real(rp)                             :: lambda_new
    real(rp)                             :: lambda_min_mine
    real(rp)                             :: lambda_max_mine

    nn    = solve % nequa
    ndof  = solve % ndofn
    nrows = solve % nequa * ndof
    ncols = solve % ncols * ndof
    nsize = max(1_ip,ncols,nrows)
    
    nullify(xx)
    nullify(yy)
    call memory_alloca(memit,'XX','solver_condition_number',xx,nsize)
    call memory_alloca(memit,'YY','solver_condition_number',yy,nsize)

    toler = 1.0e-04_rp
    !
    ! Lambda_max
    !
    xx = 1.0_rp
    call solver_parallel_vector_L2norm(solve,xx,xnorm,OPENMP=.true.)
    do ii = 1,ncols
       xx(ii) = xx(ii) / xnorm
    end do
    
    iiter  = 0
    eps    = huge(1.0_rp)
    lambda = 0.0_rp
    do while( iiter < maxit .and. eps >= toler )
       iiter = iiter + 1
       call solver_preconditioning(solve,an,xx,yy)
       !call solver_parallel_SpMV(solve,an,xx,yy,OPENMP=.true.)
       call solver_parallel_double_scalar_product(solve,yy,yy,xx,ynorm,xdoty,OPENMP=.true.)
       ynorm = max(sqrt(ynorm),zeror)
       do ii = 1,ncols
          xx(ii) = yy(ii) / ynorm
       end do
       lambda_new = xdoty 
       eps        = (lambda_new-lambda) / max(lambda_new,zeror)
       lambda     = lambda_new
       if(INOTSLAVE) print*,'max=',eps,lambda
    end do
    lambda_max_mine = lambda
    if( iiter >= maxit ) goto 10
    !
    ! Lambda_min
    !
    xx = 1.0_rp
    call solver_parallel_vector_L2norm(solve,xx,xnorm,OPENMP=.true.)
    do ii = 1,ncols
       xx(ii) = xx(ii) / xnorm
    end do
    
    iiter  = 0
    eps    = huge(1.0_rp)
    lambda = 0.0_rp
    do while( iiter < maxit .and. eps >= toler )
       iiter = iiter + 1
       call solver_preconditioning(solve,an,xx,yy)
       !call solver_parallel_SpMV(solve,an,xx,yy,OPENMP=.true.)
       do ii = 1,ncols
          yy(ii) = yy(ii) - lambda_max_mine * xx(ii)
       end do
       call solver_parallel_double_scalar_product(solve,yy,yy,xx,ynorm,xdoty,OPENMP=.true.)
       ynorm = max(sqrt(ynorm),zeror)
       do ii = 1,ncols
          xx(ii) = yy(ii) / ynorm
       end do
       lambda_new = xdoty + lambda_max_mine
       eps        = abs(lambda_new-lambda) / max(abs(lambda_new),zeror)
       lambda     = lambda_new
       if(INOTSLAVE) print*,'min=',eps,lambda
    end do
    lambda_min_mine = lambda

10  continue
    
    if( iiter >= maxit ) then
       kappa = -888.0_rp
       if( present(lambda_min) ) lambda_min = -888.0_rp
       if( present(lambda_max) ) lambda_max = -888.0_rp
    else
       kappa = abs(lambda_max_mine/lambda_min_mine)       
       if( present(lambda_min) ) lambda_min = lambda_min_mine
       if( present(lambda_max) ) lambda_max = lambda_max_mine
    end if

    call memory_deallo(memit,'XX','solver_condition_number',xx)
    call memory_deallo(memit,'YY','solver_condition_number',yy)
    
  end subroutine solver_condition_number

  !------------------------------------------------------------------------
  !> 
  !> @author  Guillaume Houzeaux
  !> @date    21/11/2017
  !> @brief   Preconditioning
  !> @details Solve system L q = A R^-1 p
  !> @}
  !------------------------------------------------------------------------

  subroutine solver_preconditioning(solve,an,xx,yy)

    type(soltyp), intent(inout)            :: solve !< Solver structure
    real(rp),     intent(in)               :: an(*) !< Matrix
    real(rp),     intent(in),     pointer  :: xx(:) !< Input
    real(rp),     intent(inout),  pointer  :: yy(:) !< Output
    real(rp),                     pointer  :: ww(:)
    integer(ip)                            :: ndof,nrows
    integer(ip)                            :: nn,nsize
    integer(ip)                            :: ncols
    real(rp)                               :: pn(2)

    nn    = solve % nequa
    ndof  = solve % ndofn
    nrows = solve % nequa * ndof
    ncols = solve % ncols * ndof
    nsize = max(1_ip,ncols,nrows)
    
    nullify(ww)
    call memory_alloca(memit,'WW','solver_preconditioning',ww,nsize)
    
    !call precon(&
    !     1_ip,ndof,nn,nrows,solve % kfl_symme,solve % kfl_preco,solve % ia,solve % ja,an,&
    !     pn,solver_invdiag,ww,xx,yy)
    call precon(&
         3_ip,ndof,nn,nrows,solve % kfl_symme,solve % kfl_preco,solve % ia,solve % ja,an,&
         pn,solver_invdiag,ww,xx,yy)
    
    call memory_deallo(memit,'WW','solver_preconditioning',ww)
    
  end subroutine solver_preconditioning

    subroutine solver_condition_number2(solve,an,kappa,lambda_min,lambda_max)

    use def_master, only : INOTSLAVE
    
    type(soltyp), intent(inout)          :: solve      !< Solver structure
    real(rp),     intent(in)             :: an(*)      !< Matrix
    real(rp),     intent(out)            :: kappa      !> Condition number
    real(rp),     intent(out),  optional :: lambda_min !> Minimum eigenvalue
    real(rp),     intent(out),  optional :: lambda_max !> Maximum eigenvalue

    integer(ip)                          :: iiter
    integer(ip)                          :: ndof,nn
    integer(ip)                          :: nrows,ncols
    integer(ip)                          :: nsize,ii
    integer(ip),  parameter              :: maxit = 100
    real(rp)                             :: toler,eps
    real(rp),     pointer                :: xx(:)
    real(rp),     pointer                :: yy(:)
    real(rp),     pointer                :: ww(:)
    real(rp),     pointer                :: rr(:)
    real(rp),     pointer                :: gg(:)
    real(rp)                             :: ynorm
    real(rp)                             :: xnorm
    real(rp)                             :: xnorm_new
    real(rp)                             :: xdoty
    real(rp)                             :: xnume,xdeno
    real(rp)                             :: lambda,alpha
    real(rp)                             :: lambda_new
    real(rp)                             :: lambda_min_mine
    real(rp)                             :: lambda_max_mine

    nn    = solve % nequa
    ndof  = solve % ndofn
    nrows = solve % nequa * ndof
    ncols = solve % ncols * ndof
    nsize = max(1_ip,ncols,nrows)
    
    nullify(xx)
    nullify(yy)
    nullify(ww)
    nullify(gg)
    nullify(rr)
    call memory_alloca(memit,'XX','solver_condition_number',xx,nsize)
    call memory_alloca(memit,'YY','solver_condition_number',yy,nsize)
    call memory_alloca(memit,'YY','solver_condition_number',ww,nsize)
    call memory_alloca(memit,'YY','solver_condition_number',gg,nsize)
    call memory_alloca(memit,'YY','solver_condition_number',rr,nsize)

    toler      = 1.0e-04_rp
    xx         = 1.0_rp
    alpha      = 1.0_rp
    lambda_new = 0.0_rp
    call solver_parallel_SpMV(solve,an,xx,ww,OPENMP=.true.)
    rr         = ww - lambda_new * xx
    call solver_preconditioning(solve,an,rr,ww) 
    
    iiter  = 0
    eps    = huge(1.0_rp)
    lambda = lambda_new
    
    do while( iiter < maxit .and. eps >= toler )
       
       iiter = iiter + 1
       
       call solver_parallel_SpMV(solve,an,xx,yy,OPENMP=.true.)
       call solver_parallel_SpMV(solve,an,ww,gg,OPENMP=.true.)
       call solver_parallel_double_scalar_product(solve,xx,xx,yy,xnorm,xdoty,OPENMP=.true.)
       lambda_new = xdoty / xnorm 

       gg         = gg - lambda_new * ww 
       call solver_parallel_double_scalar_product(solve,gg,rr,gg,xnume,xdeno,OPENMP=.true.)
       alpha      = - xnume / xdeno
       
       rr         = yy - lambda_new * xx 
       call solver_preconditioning(solve,an,rr,ww)       
       xx         = xx + alpha * ww

       eps        = abs(lambda_new-lambda) / max(abs(lambda_new),zeror)
       lambda     = lambda_new
       !if(INOTSLAVE) print*,'min=',eps,lambda,alpha
       
    end do
    lambda_min_mine = lambda

10  continue
    
    if( iiter >= maxit ) then
       kappa = -999.0_rp
       if( present(lambda_min) ) lambda_min = -999.0_rp
       if( present(lambda_max) ) lambda_max = -999.0_rp
    else
       kappa = abs(lambda_max_mine/lambda_min_mine)       
       if( present(lambda_min) ) lambda_min = lambda_min_mine
       if( present(lambda_max) ) lambda_max = lambda_max_mine
    end if

    call memory_deallo(memit,'XX','solver_condition_number',xx)
    call memory_deallo(memit,'YY','solver_condition_number',yy)
    
  end subroutine solver_condition_number2

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    13/12/2017
  !> @brief   Define groups from a field
  !> @details Define groups LGROU from a field XX
  !>
  !----------------------------------------------------------------------

  subroutine solver_define_groups_from_field(ngrou,nn,xx,lgrou,solve)
    
    integer(ip),  intent(in)                :: ngrou      !< Number of groups
    integer(ip),  intent(in)                :: nn         !< Size of the field
    real(rp),     intent(in)                :: xx(*)      !< Field
    integer(ip),  intent(inout),   pointer  :: lgrou(:)   !< Groups
    type(soltyp), intent(inout),   optional :: solve      !< Solver structure
    integer(ip)                             :: ii,igrou
    real(rp)                                :: xmini,xmaxi
    real(rp)                                :: delta,rgrou
    real(rp)                                :: dummr(2)
    real(rp),                      pointer  :: kgrou(:)

    nullify(kgrou)
    
    if( INOTMASTER ) then
       dummr(1) = -minval(xx(1:nn))
       dummr(2) =  maxval(xx(1:nn))
    end if
    call PAR_MAX(2_ip,dummr)

    xmini = -dummr(1)
    xmaxi =  dummr(2)
    delta = xmaxi-xmini 
    rgrou = real(ngrou,rp)       
    
    if( delta <= zeror ) call runend('solver_define_groups_from_field: TROUBLES!')

    if( INOTMASTER ) then
       
       call memory_deallo(memit,'LGROU','solver_define_groups_from_field',lgrou)
       call memory_alloca(memit,'LGROU','solver_define_groups_from_field',lgrou,nn)
       !
       ! Assign groups
       !
       do ii = 1,nn
          lgrou(ii) = int( ((xx(ii)-xmini-zeror) / delta )*rgrou ,ip) + 1
       end do
       !
       ! Impose fixity 
       !
       if( present(solve) ) then
          if( solve % kfl_iffix /= 0 ) then
             do ii = 1,nn
                if( solve % kfl_fixno(1,ii) > 0 ) lgrou(ii) = 0
             end do
          end if
       end if

    end if
    !
    ! Check number of groups
    !
    !allocate(kgrou(ngrou))
    !kgrou = 0
    !do ii = 1,nn
    !   igrou = lgrou(ii)
    !   if( igrou > 0 ) kgrou(igrou) = 1
    !end do
    !call PAR_SUM(kgrou)
    !ngrou = sum(kgrou)
    !deallocate(kgrou)
   
  end subroutine solver_define_groups_from_field

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    08/01/2018
  !> @brief   Warnings and error
  !> @details Detect errors in solver
  !>
  !----------------------------------------------------------------------

  subroutine solver_errors(solve)

    use def_kermod, only : kfl_ell
    use def_kermod, only : kfl_coo
    use def_kermod, only : kfl_full_rows
    use def_solver, only : SOL_SOLVER_MUMPS
    
    type(soltyp), intent(inout) :: solve      !< Solver structure
    !    
    ! Formats
    !
    if( kfl_ell == 0 .and. solve % kfl_format == SOL_ELL_FORMAT ) call runend('SOLVER_ERRORS: ACTIVATE ELL FORMAT IN .ker.dat FILE')
    if( kfl_coo == 0 .and. solve % kfl_format == SOL_COO_FORMAT ) call runend('SOLVER_ERRORS: ACTIVATE COO FORMAT IN .ker.dat FILE')
    !
    ! Full row
    !
    if( solve % kfl_full_rows /= 0 .and. kfl_full_rows == 0 ) &
         call runend('SOLVER_ERRORS: ACTIVATE THE FULL_ROWS OPTION IN ker.dat FILE')
    !
    ! Full row
    !
    if( solve % kfl_full_rows /= 0 .and. solve % kfl_format /= SOL_CSR_FORMAT ) &
         call runend('SOLVER_ERRORS: FULL_ROWS OPTION IS ONLY AVAILABLE WITH CSR MATRIX FORMAT')
    !
    ! Full row graph must be activated for 
    !
    if( solve % kfl_algso == SOL_SOLVER_MUMPS ) then
       if( kfl_full_rows == 0 ) &
            call runend('SOLVER_ERRORS: ACTIVATE THE FULL_ROWS OPTION IN ker.dat FILE TO USE MUMPS')
#ifndef MUMPS
       call runend('SOLVER_ERRORS: COMPILE WITH OPTION -DMUMPS TO USE MUMPS')
#endif
    end if

  end subroutine solver_errors

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-03
  !> @brief   SpMV
  !> @details Parallel SpMV
  !>          Multiply a matrix by a vector YY = A XX
  !>          where A is in CSR, COO or ELL format. This is the parallel
  !>          version, including blocking and non-blocking send-receive
  !>          MPI functions.
  !>
  !>          INPUT
  !>             NBNODES .... Number of equations
  !>             NBVAR ...... Number of variables
  !>             AN ......... Matrix
  !>             JA ......... List of elements
  !>             IA ......... Pointer to list of elements
  !>             XX ......... Vector
  !>          OUTPUT
  !>             YY ......... result vector
  !>
  !> @}
  !------------------------------------------------------------------------

  subroutine solver_SpMV(solve,an,xx,yy,MPI)

    type(soltyp), intent(inout)         :: solve  !< Solver type
    real(rp),     intent(in)            :: an(*)  !< Matrix
    real(rp),     intent(inout)         :: xx(*)  !< Input vector
    real(rp),     intent(out)           :: yy(*) !< Output vector
    logical(lg),  intent(in),  optional :: MPI   !< If MPI should be used or not

#ifdef EXTRAE
    call extrae_eventandcounters(900,int(5,8)) ! Matrix-vector product
#endif

    if( solve % omp_interface == 0 ) then
       call solver_parallel_SpMV(solve,an,xx,yy,OPENMP=.true.,TIMING=.true.,OPENMP_INTERFACE=.false.,MPI=MPI)
    else
       call solver_parallel_SpMV(solve,an,xx,yy,OPENMP=.true.,TIMING=.true.,MPI=MPI)
    end if

#ifdef EXTRAE
    call extrae_eventandcounters(900,int(0,8)) ! Matrix-vector product
#endif

  end subroutine solver_SpMV

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-04
  !> @brief   Dirichlet conditions
  !> @details Impose Dirichlet conditions on an unknown
  !> 
  !-----------------------------------------------------------------------

  subroutine solver_impose_dirichlet_condition(solve,xx)

    type(soltyp), intent(inout) :: solve    !< Solver type
    real(rp),     intent(out)   :: xx(*)    !< Unknown
    integer(ip)                 :: ii,kk
    integer(ip)                 :: ndofn
    
    ndofn = solve % ndofn
    if( solve % kfl_iffix == 1 .and. associated(solve % bvess) ) then
       !
       ! Dirichlet condition prescribed in BVESS
       !
       do ii = 1,solve % nequa
          do kk = 1,ndofn
             if( solve % kfl_fixno(kk,ii) > 0 ) xx((ii-1)*ndofn+kk) = solve % bvess(kk,ii)
          end do
       end do

    else if( solve % kfl_iffix == 2 .or. ( solve % kfl_iffix == 1 .and. .not. associated(solve % bvess) ) ) then
       !
       ! Dirichlet condition = 0
       !
       do ii = 1,solve % nequa
          do kk = 1,ndofn
             if( solve % kfl_fixno(kk,ii) > 0 ) xx((ii-1)*ndofn+kk) = 0.0_rp
          end do
       end do
       
    end if
    
  end subroutine solver_impose_dirichlet_condition

end module mod_solver
!> @}
!-----------------------------------------------------------------------
