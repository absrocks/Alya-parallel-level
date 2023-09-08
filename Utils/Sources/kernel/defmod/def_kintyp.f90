!-----------------------------------------------------------------------
!> @defgroup Kinds_and_types
!> Kinds ands types of Alya
!> @{
!> @file    def_kintyp.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   Definition of kinds and types.
!> @details Definition of kinds and types.
!>          "The range of the default integers is not specified in the language
!>          but on a computer with a word size of n bits, is often from
!>          -2^{n-1} to +2^{n-1}-1. Thus on a 32-bit computer the range is
!>          often -2.14*10^9 to +2.14*10^9."
!>          M. Metclaf and J. Reid, FORTRAN 90/95 explained, 2nd edition.
!>
!>          Defaults are:
!>          Integers: 4-bytes
!>          Reals:    8-bytes
!>
!-----------------------------------------------------------------------

module def_kintyp

#ifdef MAPHYS
  use DMPH_maphys_mod
#endif
#ifdef MUMPS
  INCLUDE 'dmumps_struc.h'
#endif

  !----------------------------------------------------------------------
  !
  ! Symbolc names for integers, reals and logicals
  !
  !----------------------------------------------------------------------
  !
  ! Symbolic names for integers
  !
#ifdef I8
  integer, parameter  :: ip = 8             ! 8-byte integer
#else
  integer, parameter  :: ip = 4             ! 4-byte integer
#endif
  !
  ! Symbolic names for reals
  !
#ifdef R4
  integer, parameter  :: rp = 4             ! Simple precision
#elif R16
  integer, parameter  :: rp = 16            ! Very high precision - Beware does not work with MPI
#else
  integer, parameter  :: rp = 8             ! Double precision
#endif
  !
  ! Symbolic name for kind type of default logical
  !
  integer, parameter  :: lg = kind(.true.)

  !----------------------------------------------------------------------
  !
  ! General types
  !
  !----------------------------------------------------------------------

  type i1p
     integer(ip), pointer :: l(:)
  end type i1p
  type i2p
     integer(ip), pointer :: l(:,:)
  end type i2p
  type i3p
     integer(ip), pointer :: l(:,:,:)
  end type i3p
  type r1p
     real(rp),    pointer :: a(:)
  end type r1p
  type r2p
     real(rp),    pointer :: a(:,:)
  end type r2p
  type r3p
     real(rp),    pointer :: a(:,:,:)
  end type r3p
  type r4p
     real(rp),    pointer :: a(:,:,:,:)
  end type r4p
  type i1pp
     integer(ip)          :: n
     integer(ip), pointer :: l(:)
  end type i1pp
  type i1pi1p
     type(i1p),   pointer :: l(:)
  end type i1pi1p
  type spmat
     integer(ip)          :: ndof
     integer(ip)          :: nrows
     integer(ip)          :: ncols
     integer(ip), pointer :: iA(:)
     integer(ip), pointer :: jA(:)
     real(rp),    pointer :: vA(:,:,:)
  end type spmat
  !----------------------------------------------------------------------
  !
  ! Boundary codes
  !
  !----------------------------------------------------------------------

  type bccodes
     integer(ip)          :: varna
     real(rp)             :: xcent(3)
     integer(ip)          :: ldofr(5)
     integer(ip)          :: funty
     real(rp)             :: param(6)
  end type bccodes
  !
  ! Node codes
  !
  type bc_nodes1
     integer(ip),     pointer :: lcode(:)
     integer(ip)              :: kfl_fixno
     integer(ip)              :: kfl_value
     real(rp),        pointer :: bvess(:)
     integer(ip)              :: kfl_funno ! Old time functions
     integer(ip)              :: kfl_fixrs
     character(5)             :: tag
     character(5)             :: fname
  end type bc_nodes1
  type bc_nodes
     integer(ip)              :: kfl_ibopo
     integer(ip)              :: ndofn
     integer(ip)              :: ncode
     type(bc_nodes1), pointer :: l(:)
  end type bc_nodes
  !
  ! Boundary codes
  !
  type bc_bound1
     integer(ip)              :: lcode
     integer(ip)              :: kfl_fixbo
     integer(ip)              :: kfl_value
     integer(ip)              :: kfl_funbo
     character(5)             :: tag
     real(rp),        pointer :: bvnat(:)
     character(5)             :: fname
  end type bc_bound1
  type bc_bound
     integer(ip)              :: ndofn
     integer(ip)              :: ncode
     type(bc_bound1), pointer :: l(:)
  end type bc_bound
  !
  ! Space time function
  !
  type typ_space_time_function
     integer(ip)             :: ndime         ! Number of dimensions
     integer(ip)             :: nexpr         ! Size of the expression
     integer(ip)             :: numfield      ! Field number, when the function is a field
     character(5)            :: name          ! Name
     character(400), pointer :: expression(:) ! Expression to be parsed
  end type typ_space_time_function
  !
  ! Time function
  !
  type typ_time_function
     integer(ip)             :: kfl_type      ! Function type
     integer(ip)             :: npara         ! Number of parameters
     real(rp),       pointer :: parameters(:) ! Parameters
     character(5)            :: name          ! Name
  end type typ_time_function

  !----------------------------------------------------------------------
  !
  ! Element data base type
  !
  !----------------------------------------------------------------------

  type elm
     ! User integration rule
     integer(ip)          :: pgaus            ! Number of Gauss points
     real(rp),    pointer :: shape(:,:)       ! pnode,pgaus
     real(rp),    pointer :: deriv(:,:,:)     ! ndime,pnode,pgaus
     real(rp),    pointer :: heslo(:,:,:)     ! ntens,pnode,pgaus
     real(rp),    pointer :: weigp(:)
     real(rp),    pointer :: shaga(:,:)
     real(rp),    pointer :: posgp(:,:)
     ! Bubble
     real(rp),    pointer :: shape_bub(:)     ! pgaus
     real(rp),    pointer :: deriv_bub(:,:)   ! ndime,pgaus
     real(rp),    pointer :: heslo_bub(:,:)   ! ntens,pgaus
     ! Center of gravity
     real(rp),    pointer :: shacg(:)
     real(rp),    pointer :: dercg(:,:)
     real(rp),    pointer :: hescg(:,:)
     real(rp)             :: weicg
     ! Close rule
     real(rp),    pointer :: shapc(:,:)
     real(rp),    pointer :: deric(:,:,:)
     real(rp),    pointer :: heslc(:,:,:)
     real(rp),    pointer :: weigc(:)
     ! IB integration rule
     real(rp),    pointer :: shaib(:,:)
     real(rp),    pointer :: derib(:,:,:)
     real(rp),    pointer :: weiib(:)
  end type elm

  type elmgp
     real(rp),    pointer :: gpvol(:)
     real(rp),    pointer :: gpcar(:,:,:)
     real(rp),    pointer :: gphes(:,:,:)
     real(rp),    pointer :: hleng(:)
     real(rp),    pointer :: tragl(:,:)
  end type elmgp

  !----------------------------------------------------------------------
  !
  ! Communication arrays
  !
  !----------------------------------------------------------------------

  type comm_bound_matrix
     type(i1p),   pointer :: ja(:)
     integer(ip), pointer :: nzdom_ii(:)
     integer(ip)          :: nzdom
  end type comm_bound_matrix
  type comm_data_par
     integer(ip)                      :: nneig
     integer(ip), pointer             :: neights(:)
     ! Interior, own and other boundary nodes
     integer(ip)                      :: npoi1
     integer(ip)                      :: npoi2
     integer(ip)                      :: npoi3
     ! Interior, own and other boundary edges
     integer(ip)                      :: nedg1
     integer(ip)                      :: nedg2
     integer(ip)                      :: nedg3
     ! Interface node communication
     integer(ip), pointer             :: bound_size(:)
     integer(ip), pointer             :: bound_perm(:)
     integer(ip), pointer             :: bound_scal(:)
     integer(ip)                      :: bound_dim
     integer(ip), pointer             :: bound_multiplicity(:)
     integer(ip), pointer             :: bound_owner_rank(:)
     integer(ip), pointer             :: node_number_in_owner(:)
     type(comm_bound_matrix), pointer :: bound_matrix(:)
     type(comm_bound_matrix), pointer :: bound_mat_halo_send(:)
     type(comm_bound_matrix), pointer :: bound_mat_halo_recv(:)
     ! Interface edge communication
     integer(ip), pointer             :: bedge_size(:)
     integer(ip), pointer             :: bedge_perm(:)
     integer(ip), pointer             :: bedge_adja(:)
     integer(ip), pointer             :: bedge_scal(:)
     integer(ip)                      :: bedge_dim
     integer(ip), pointer             :: bedge_multiplicity(:)
     integer(ip), pointer             :: bedge_owner_rank(:)
     ! Non-symmetric send receive
     integer(ip), pointer             :: lsend_size(:)
     integer(ip), pointer             :: lsend_perm(:)
     integer(ip)                      :: lsend_dim
     integer(ip), pointer             :: lrecv_size(:)
     integer(ip), pointer             :: lrecv_perm(:)
     integer(ip)                      :: lrecv_dim
     integer(ip), pointer             :: lscat_perm(:)
     integer(ip)                      :: lscat_dim
     integer(ip)                      :: matrix_nzdom
     integer(ip), pointer             :: matrix_ia(:)
     integer(ip), pointer             :: matrix_ja(:)
     real(rp),    pointer             :: matrix_aa(:)
     ! Face communication
     integer(ip), pointer             :: bface_size(:)
     integer(ip), pointer             :: bface_perm(:)
     integer(ip)                      :: bface_dim
     ! Full row matrix communicator
     integer(ip)                      :: full_row_send_nneig
     integer(ip), pointer             :: full_row_send_neights(:)
     integer(ip), pointer             :: full_row_send_size(:)
     integer(ip), pointer             :: full_row_send_perm(:)
     integer(ip)                      :: full_row_send_dim
     integer(ip)                      :: full_row_recv_nneig
     integer(ip), pointer             :: full_row_recv_neights(:)
     integer(ip), pointer             :: full_row_recv_size(:)
     integer(ip), pointer             :: full_row_recv_perm(:)
     integer(ip)                      :: full_row_recv_dim
     ! Ghost node communication
     integer(ip), pointer             :: ghost_send_node_size(:)
     integer(ip), pointer             :: ghost_send_node_perm(:)
     integer(ip)                      :: ghost_send_node_dim
     integer(ip), pointer             :: ghost_recv_node_size(:)
     integer(ip), pointer             :: ghost_recv_node_perm(:)
     integer(ip)                      :: ghost_recv_node_dim
     ! Ghost element communication
     integer(ip), pointer             :: ghost_send_elem_size(:)
     integer(ip), pointer             :: ghost_send_elem_perm(:)
     integer(ip)                      :: ghost_send_elem_dim
     integer(ip), pointer             :: ghost_recv_elem_size(:)
     integer(ip), pointer             :: ghost_recv_elem_perm(:)
     integer(ip)                      :: ghost_recv_elem_dim
     ! Ghost boundary communication
     integer(ip), pointer             :: ghost_send_boun_size(:)
     integer(ip), pointer             :: ghost_send_boun_perm(:)
     integer(ip)                      :: ghost_send_boun_dim
     integer(ip), pointer             :: ghost_recv_boun_size(:)
     integer(ip), pointer             :: ghost_recv_boun_perm(:)
     integer(ip)                      :: ghost_recv_boun_dim
     ! Communicator
     integer(ip)                      :: PAR_COMM_WORLD
     integer(4)                       :: RANK4
  end type comm_data_par

  !----------------------------------------------------------------------
  !
  ! Solver
  !
  !----------------------------------------------------------------------

#ifndef MAPHYS
  type DMPH_maphys_t
     integer(ip)          :: n
     integer(ip)          :: nnz
     integer(ip)          :: sym
     integer(ip)          :: job
     integer(4)           :: comm
     real(rp),    pointer :: values(:)
     real(rp),    pointer :: rhs(:)
     real(rp),    pointer :: sol(:)
     real(rp),    pointer :: RCNTL(:)
     real(rp),    pointer :: RINFOG(:)
     real(rp),    pointer :: RINFO(:)
     integer(ip), pointer :: IINFOG(:)
     integer(ip), pointer :: ICNTL(:)
     integer(ip), pointer :: rows(:)
     integer(ip), pointer :: cols(:)
  end type DMPH_maphys_t
#endif

  type lpoin_block_typ
     type(block_matrix_typ), pointer :: block2_num(:,:)
     type(block_rhs_typ),    pointer :: block1_num(:)
  end type lpoin_block_typ
  type block_matrix_typ
     real(rp), pointer :: matrix(:,:,:)
  end type block_matrix_typ
  type block_rhs_typ
     real(rp),    pointer :: reaction(:)
     real(rp),    pointer :: rhs(:)
     real(rp),    pointer :: bvess(:,:)
     integer(ip), pointer :: kfl_fixno(:,:)     ! Fixity array
     real(rp),    pointer :: bvnat(:,:)         ! Neumann value
  end type block_rhs_typ
  type block_typ
     real(rp),    pointer :: bvess(:,:)
     integer(ip), pointer :: kfl_fixno(:,:)     ! Fixity array
     real(rp),    pointer :: bvnat(:,:)         ! Neumann value
  end type block_typ

#ifdef PASTIX
#include "pastix_fortran.h"
#else
!  #define PASTIX_INT_KIND    4
!  #define pastix_int_t       INTEGER(kind=4)
!  #define pastix_uint_t      unsigned INTEGER(kind=4)
!  #define pastix_data_ptr_t  INTEGER(kind=8)
!  #define MPI_PASTIX_INT     MPI_INTEGER4
!  #define pastix_float_t     REAL(kind=8)
!  #define MPI_PASTIX_FLOAT   MPI_REAL8
#endif

  type direct_solver_typ
     ! Input data
     integer(ip)          :: kfl_solver          ! Solver type
     integer(ip)          :: kfl_symmetric       ! Symmetric matrix
     integer(ip)          :: nn                  ! # nodes
     integer(ip)          :: nn_own              ! # own nodes
     integer(ip)          :: ndof                ! # number dof per node
     integer(ip)          :: nrhs                ! # RHS
     integer(ip)          :: nz                  ! CSR or CSC: # non-zero edges
     integer(ip)          :: nz_own              ! CSR or CSC: # non-zero edges
     integer(ip)          :: kfl_paral           ! Should be run in parallel
     integer(ip), pointer :: ia(:)               ! CSR or CSC: Matrix graph
     integer(ip), pointer :: ja(:)               ! CSR or CSC: Matrix graph
     integer(ip)          :: nskyl               ! Skyline: matrix size
     integer(ip), pointer :: idiag(:)            ! Skyline: diagonal
     integer(ip), pointer :: iskyl(:)            ! Skyline: pointer
     ! Output
     integer(8)           :: memor(2)            ! Memory counter
     real(rp)             :: cputi(3)            ! CPU time
     integer(ip)          :: num_initializations ! Number of initializations
     integer(ip)          :: num_solutions       ! Number of solves
     integer(ip)          :: num_factorizations  ! Number of factorizations
     character(30)        :: name                ! Name of the solver
     real(rp)             :: fillin              ! Fill-in
     ! Solver internal arrays
     real(rp),    pointer :: aa(:,:,:)           ! Matrix
     integer(ip), pointer :: ia_new(:)           ! Alya renumbered graph
     integer(ip), pointer :: ja_new(:)           ! Alya renumbered graph
     real(rp),    pointer :: aa_skyl(:)          ! Skyline factorized matrix
     integer(ip), pointer :: IL(:)
     integer(ip), pointer :: JL(:)
     real(rp),    pointer :: LN(:)
     integer(ip), pointer :: IU(:)
     integer(ip), pointer :: JU(:)
     real(rp),    pointer :: UN(:)
     integer(ip), pointer :: permr(:)            ! For Pastix, pastix_int_t
     integer(ip), pointer :: invpr(:)            ! For Pastix, pastix_int_t
     ! Pastix internal arrays
     integer(ip), pointer :: iparm(:)            ! For Pastix, pastix_int_t, WSMP and PWSMP
     real(rp),    pointer :: dparm(:)            ! For Pastix, pastix_float_t,WSMP and PWSMP
     integer(ip), pointer :: bindtab(:)
     integer(8)           :: pastix_data         ! For Pastix, pastix_data_ptr_t
     integer(4)           :: pastix_comm         ! For Pastix, pastix_mpi_int
     integer(ip)          :: ldb
     integer(ip)          :: num_threads
     real(rp),    pointer :: avals(:)
     ! MUMPS
     integer(ip), pointer :: gatsca_mumps(:)     ! Gather MUMP
#ifdef MUMPS
     TYPE (DMUMPS_STRUC) mumps_par
#endif
  end type direct_solver_typ

  type soltyp
     !
     ! Input
     !
     integer(ip)          :: kfl_algso          ! Solver tag
     integer(ip)          :: ndofn              ! # dof per nodes
     integer(ip)          :: kfl_recov          ! Recover original residual after slave exchange
     integer(ip)          :: kfl_ortho          ! Orthogonolization GMRES (0=classical Gram-Schmidt,1=modified)
     integer(ip)          :: kfl_limit          ! Algebraic limiter
     integer(ip)          :: kfl_adres          ! Adaptive residual (0=no,1=yes)
     integer(ip)          :: nrhss              ! # of simultaneous RHS
     integer(ip)          :: miter              ! Maximum number of iterations
     integer(ip)          :: nkryd              ! Krylov dimension
     integer(ip)          :: kfl_coarse         ! If a coarse solver is to be used
     real(rp)             :: solco              ! Solver tolerance
     real(rp)             :: adres              ! Adaptive residual tolerance (0=no,1=yes)
     real(rp)             :: solmi              ! Minimum solver tolerance (if adaptive tolerance)
     integer(ip)          :: kfl_cvgso          ! Convergence output
     integer(ip)          :: lun_cvgso          ! Convergence file unit
     integer(ip)          :: kfl_solve          ! Solver info
     integer(ip)          :: lun_solve          ! Solver info file unit
     integer(ip)          :: lun_exsol          ! External solver info file unit
     integer(ip)          :: kfl_penal          ! Penalization: method
     real(rp)             :: penal              ! Penalization: parameter
     integer(ip)          :: kfl_defpr          ! Smoother preconditioner
     integer(ip)          :: itpre              ! Jacobi iterations
     integer(ip)          :: kfl_preco          ! Preconditioner tag
     integer(ip)          :: kfl_leftr          ! Left(0) or right(1)
     integer(ip)          :: kfl_marke          ! Output of matrix in market format
     integer(ip)          :: kfl_force          ! Force continuity of solution across interface after solver
     integer(ip)          :: kfl_format         ! CSR (1), COO (2), ELL (3)
     integer(ip)          :: kfl_where          ! Unknown is on node, edge, element center
     integer(ip)          :: kfl_clean_precond  ! If preconditioner/coarse solver should be always recomputed
     integer(ip)          :: kfl_normalization  ! Residual normalization strategy
     integer(ip)          :: kfl_block_ras      ! Full or block RAS
     integer(ip)          :: num_multiple_solves! Number of multiple solves using the same matrix
     integer(ip)          :: kfl_full_rows      ! Partial or full row matrix
     integer(ip)          :: kfl_save_krylov    ! If Krylov subspace should be saved
     integer(ip)          :: kfl_kappa          ! If condition number should be computed
     integer(ip)          :: kfl_roe_correction ! Round off error corrections should be used
     real(rp)             :: normalization      ! Residual normalization value
     character(50)        :: conf_file          ! Configuration file
     real(rp)             :: bnorm_min          ! Minimum norm of RHS to consider null system
     real(rp)             :: threshold          ! Threshold used for direct solver base preconditioners
     integer(ip)          :: omp_schedule       ! OMP schedule
     integer(ip)          :: omp_chunk_size     ! OMP chunk size
     integer(ip)          :: omp_interface      ! OMP for SpMV on interface
     integer(ip)          :: kfl_iffix          ! If solver takes care of imposing Dirichlet b.c.
     integer(ip)          :: kfl_bvnat          ! If natural b.c. are imposed in solver
     integer(ip)          :: kfl_dirichlet      ! Way dirichlet b.c. are imposed
     !
     ! Output
     !
     integer(ip)          :: iters              ! # iterations
     integer(ip)          :: itsol(3)           ! Solver statistics
     integer(ip)          :: nsolv              ! # solves
     real(rp)             :: resin              ! Initial preconditioned residual
     real(rp)             :: resfi              ! Final preconditioned residual
     real(rp)             :: resi2              ! Initial non-preconditioned residual = ||b-Ax||/||b||
     real(rp)             :: resf2              ! Final non-preconditioned residual = ||b-Ax||/||b||
     real(rp)             :: cputi(10)          ! CPU time
     integer(ip)          :: num_spmv           ! Number of SpMV
     integer(ip)          :: num_dot            ! Number of dot products
     real(rp)             :: cpu_spmv(9)        ! CPU time of SMVP
     real(rp)             :: cpu_dot(9)         ! CPU time of dot product
     real(rp)             :: cpu_schur(10)      ! Schur complement solver timing
     real(rp)             :: xorth              ! Orthogonality check = (p^i-1,p^i)
     integer(ip)          :: kfl_exres          ! Exact residual output ||b-Ax||/||b||
     real(rp)             :: xdiag              ! Coefficient for diagonal solver
     integer(ip)          :: num_solves         ! Number of solves for multiple solves
     real(rp)             :: bnorm              ! RHS L2 norm
     real(rp)             :: kappa              ! Condition number of the matrix
     real(rp)             :: lambda_min         ! Minimum eigenvalue
     real(rp)             :: lambda_max         ! maximum eigenvalue
     !
     ! Derived parameters
     !
     character(50)        :: wprob              ! Problem name
     character(50)        :: wsolv              ! Solver name
     character(50)        :: wprec              ! Preconditioner name
     integer(ip)          :: nzmat              ! System size
     integer(ip)          :: nzmat_own          ! System size for full row matrix
     integer(ip)          :: nzrhs              ! Size RHS = ndofn * nequa
     integer(ip)          :: nequa              ! # nodes to solve
     integer(ip)          :: nequa_own          ! # nodes to solve (own nodes for scalar product)
     integer(ip)          :: nequa_halo         ! # nodes to solve (own nodes for scalar product)
     integer(ip)          :: nunkn              ! # unknowns to solve nequa*ndofn
     integer(ip)          :: ncols              ! # columns (non square matrix)
     integer(ip)          :: ndof2              ! # dof^2
     integer(ip)          :: nequ1              ! Interior unknowns
     integer(ip)          :: nequ2              ! Start own boundary
     integer(ip)          :: nequ3              ! End own boundary
     integer(ip)          :: kfl_symme          ! Matrix symmetric assembly (0=no,1=yes)
     integer(ip)          :: kfl_symeq          ! Equation symmetry (0=no,1=yes)
     integer(ip)          :: kfl_cmplx          ! Real or complex solver JELENA
     integer(ip)          :: kfl_assem          ! Matrix has been assembled
     integer(ip)          :: kfl_schur          ! Schur solver (0=no)
     integer(ip)          :: kfl_schum          ! If matrices come from a Schur complement
     integer(ip)          :: kfl_version        ! Version of solver (new one incldues pre and post)
     integer(ip)          :: heade              ! Track if header has already been written in file
     character(150)       :: fil_cvgso          ! Convergence file
     character(150)       :: fil_solve          ! Solver file
     integer(ip)          :: nzpre              ! Size of preconditioner matrix
     integer(ip)          :: kfl_update_precond ! If preconditionner should be updated
     !
     ! Mask for dot product
     !
     integer(ip)          :: kfl_mask           ! If mask
     real(rp),    pointer :: mask(:)            ! The mask
     !
     ! Matrix graph
     !
     integer(ip)          :: nnz                ! Size of graph
     integer(ip)          :: nnz_ell            ! Size of ELL graph
     integer(ip), pointer :: ia(:)              ! CSR format: IA
     integer(ip), pointer :: ja(:)              ! CSR format: JA
     integer(ip), pointer :: ia_full(:)         ! CSR format: IA for full row format
     integer(ip), pointer :: ia_full_end(:)     ! CSR format: IA for full row format, end until own nodes
     integer(ip), pointer :: ia_full_ini(:)     ! CSR format: IA for full row format, start after own nodes
     integer(ip), pointer :: ja_full(:)         ! CSR format: JA for full row format
     integer(ip), pointer :: rows(:)            ! COO format: rows
     integer(ip), pointer :: cols(:)            ! COO format: columns
     integer(ip), pointer :: cols_ell(:,:)      ! ELL format: columns
     !
     ! Block Gauss-Seidel
     !
     integer(ip)          :: kfl_blogs          ! Block Gauss-Seidel treatment
     integer(ip)          :: nblok              ! Number of blocks
     integer(ip)          :: ndofn_per_block    ! Number of degress of freedom per block
     integer(ip), pointer :: lperm_block(:,:)   ! Permutation arrays BGS
     integer(ip), pointer :: linvp_block(:,:)   ! Inverse permutation arrays BGS
     !
     ! Fixity and boundary conditions
     !
     real(rp),    pointer :: bvess(:,:)         ! Dirichlet value
     integer(ip), pointer :: kfl_fixno(:,:)     ! Fixity array
     real(rp),    pointer :: bvnat(:,:)         ! Neumann value
     !
     ! Direct solver
     !
     type(direct_solver_typ)          :: direct_solver              ! Direct solver
     type(direct_solver_typ), pointer :: direct_solver_RAS(:)       ! Direct solver for RAS
     type(direct_solver_typ)          :: direct_solver_coarse       ! Direct solver for coarse
     type(direct_solver_typ)          :: direct_solver_block_LU     ! Direct solver for block LU
     type(direct_solver_typ)          :: direct_solver_AMG          ! Direct solver for AMG
     type(direct_solver_typ)          :: direct_solver_Deflation    ! Direct solver for deflation
     !
     ! Deflated CG: deflcg.f90
     !
     integer(ip)          :: ngrou              ! Deflated CG: # groups
     integer(ip)          :: nskyl              ! Deflated CG: Skyline A' matrix size
     integer(ip)          :: ifbop              ! Deflated CG: If boundary condition is imposed on boundary nodes
     integer(ip)          :: kfl_gathe          ! Deflated CG: All reduce / or gather
     integer(ip)          :: kfl_defas          ! Deflated CG: Assembly of deflated
     integer(ip)          :: kfl_defso          ! Deflated CG: Solver of deflated
     integer(ip), pointer :: limpo(:)           ! Deflated CG: List of imposed nodes
     integer(ip), pointer :: lgrou(:)           ! Deflated CG: Group list
     integer(ip), pointer :: iskyl(:)           ! Deflated CG: Skyline index
     integer(ip), pointer :: idiag(:)           ! Deflated CG: Pointer to skyline diagonal
     integer(ip), pointer :: iagro(:)           ! Deflated CG: Sparse graph
     integer(ip), pointer :: jagro(:)           ! Deflated CG: Sparse graph
     integer(ip)          :: nzgro              ! Deflated CG: Sparse graph
     integer(ip)          :: icoml              ! Deflated CG: Communication level (for parallelization)
     integer(ip)          :: nbig               ! Deflated CG: All gather strategy
     integer(ip)          :: nsmall             ! Deflated CG: All gather strategy
     integer(ip), pointer :: lcoun(:)           ! Deflated CG: All gather strategy
     integer(ip), pointer :: lbig(:)            ! Deflated CG: All gather strategy
     integer(ip), pointer :: displ(:)           ! Deflated CG: All gather strategy
     integer(4),  pointer :: disp4(:)           ! Deflated CG: All gather strategy
     integer(4),  pointer :: lcou4(:)           ! Deflated CG: All gather strategy
     real(rp),    pointer :: xsmall(:)          ! Deflated CG: All gather strategy
     real(rp),    pointer :: xbig(:)            ! Deflated CG: All gather strategy
     !
     ! Renumbered Gauss-Seidel
     !
     integer(ip)          :: kfl_renumbered_gs  ! Renumbered Gauss-Seidel permutation
     integer(ip), pointer :: permr_gs(:)        ! Renumbered Gauss-Seidel permutation
     integer(ip)          :: ngrou_gs           ! Renumbered Gauss-Seidel group number
     integer(ip), pointer :: invpr_gs(:)        ! Renumbered Gauss-Seidel inverse permutation
     integer(ip), pointer :: lgrou_gs(:)        ! Renumbered Gauss-Seidel groups
     real(rp),    pointer :: vecto_gs(:,:)      ! Renumbered Gauss-Seidel vetor
     integer(ip), pointer :: idiag1(:,:)        ! Stores the position of the diagonal terms in Streamwise Bidiagonal precon.
     real(rp),    pointer :: adiag1(:,:,:)      ! Stores the subdiagonal matrix coefficents in Streamwise Bidiagonal precon.
     real(rp)             :: angle_stream       ! Minimum angle to determine next node renumbered in Streamwise direction
     !
     ! Linelet
     !
     integer(ip)          :: npntr              ! Linelet: # non-zero in trima
     integer(ip)          :: nlpntr             ! Linelet: # points in linelet
     integer(ip)          :: nline              ! Linelet: # linelets
     integer(ip)          :: nlin1              ! Linelet: # linelets crossed by boundary
     integer(ip)          :: npoin              ! Linelet: % of nodes in linelets
     integer(ip)          :: kfl_factl          ! Linelet: Factorization flag
     integer(ip)          :: kfl_linty          ! Linelet: Linlet type
     integer(ip), pointer :: lpntr(:)           ! Linelet: tridiag. to original sym. matrix
     integer(ip), pointer :: lrenu(:)           ! Linelet: point position in tridiag.
     integer(ip), pointer :: lrenup(:)          ! Linelet: inverse permutation
     integer(ip), pointer :: limli(:)           ! Linelet: list of imposed node
     integer(ip), pointer :: lline(:)           ! Linelet: Pointer for each linelet
     real(rp),    pointer :: trima(:)           ! Linelet: Tridiagonal matrix
     real(rp)             :: toler              ! Linelet: Tolerance aspect ratio
     !
     ! Schur solver
     !
     integer(ip)          :: kfl_scpre          ! Schur preconditioner solver
     integer(ip)          :: kfl_scaii          ! Schur matrix solver
     integer(ip)          :: poaii              ! Pointer to Aii
     integer(ip)          :: poaib              ! Pointer to Aib
     integer(ip)          :: poabi              ! Pointer to Abi
     integer(ip)          :: poabb              ! Pointer to Abb
     !
     ! Block character of matrix, used to compute Reaction residuals
     ! Save all blocks of matrix
     ! e.g.: block(1,2) % a(:,:,:) = matrix of block 1,2 in CSR format
     !       block_number = 3
     !       block_dimensions(1:3) = 1,ndime,1
     !
     ! +-----+--------+-----+
     ! |     |        |     |
     ! | 1,1 |  1,2   | 1,3 |
     ! +-----+--------+-----+
     ! |     |        |     |
     ! | 2,1 |  2,2   | 2,3 |
     ! |     |        |     |
     ! +-----+--------+-----+
     ! |     |        |     |
     ! | 3,1 |  3,2   | 3,3 |
     ! +-----+--------+-----+
     !
     ! lpoin_block(1:npoin) % block2_num(i,j) % matrix(:,:,:)
     ! lpoin_block(1:npoin) % block2_num(i,j) % rhsi(:)
     ! lpoin_block(1:npoin) % block1_num(i)   % reaction(:)
     !
     integer(ip)                     :: kfl_react              ! If reaction should be computed
     type(lpoin_block_typ), pointer  :: lpoin_block(:)
     logical(lg),           pointer  :: lpoin_reaction(:)
     real(rp),              pointer  :: reaction(:,:)
     integer(ip)                     :: num_blocks
     integer(ip)                     :: block_num
     integer(ip)                     :: block_dimensions(10)
     integer(ip)                     :: block_pointers(10)
     type(block_typ)                 :: block_array(10)
     !
     ! Main matrix comes from a Schur complement
     ! A = A1 + A2*A3^-1*A4
     !
     integer(ip)          :: ndofn_A3
     real(rp),    pointer :: A1(:)
     real(rp),    pointer :: A2(:)
     real(rp),    pointer :: A3(:)
     real(rp),    pointer :: A4(:)
     real(rp),    pointer :: invA3(:)
     !
     ! Krylov subspace
     !
     real(rp),    pointer :: krylov_dir(:,:)
     real(rp),    pointer :: krylov_dot_products(:)
     integer(ip)          :: krylov_size
     !
     ! External solvers
     !
     type(DMPH_maphys_t)  :: mphs                 ! MAPHYS or dummy type

  end type soltyp

  !----------------------------------------------------------------------
  !
  ! Eigenvalue solver
  !
  !----------------------------------------------------------------------

  type eigtyp
     character(50)        :: wprob              ! Names
     character(50)        :: wsolv
     character(50)        :: wprec
     integer(ip)          :: ndofn              ! D.o.F
     integer(ip)          :: ndof2              ! D.o.F*D.o.F
     integer(ip)          :: neiva              ! number of eigen values requiered
     integer(ip)          :: neige              ! Size of eigenvector
     integer(ip)          :: nzmbt              ! RHS eigen matrix
     integer(ip)          :: nzmat              ! LHS matrix
     integer(ip)          :: nzrhs              ! RHS
     integer(ip)          :: kfl_algso          ! Solver parameter
     integer(ip)          :: kfl_massm          ! Mass matrix type
     integer(ip)          :: miter              ! Solver number of iterations
     integer(ip)          :: kfl_facto          ! Factorization
     integer(ip)          :: itsol(3)           ! Output: Solver statistics
     integer(ip)          :: nsolv              ! Output: # solves
     integer(ip)          :: lun_cvgei
     integer(ip)          :: lun_solei
     integer(ip)          :: kfl_preco          ! Preconditioner
     real(rp)             :: solco              ! Tolerance
     real(rp)             :: shift              ! shift
     integer(ip)          :: diter              ! iter en eigdir
     integer(ip)          :: error              ! error en eigdir
     integer(ip), pointer :: lpdof(:)
  end type eigtyp

  !----------------------------------------------------------------------
  !
  ! atomos y especies pointers
  !
  !----------------------------------------------------------------------

  type atomo
     character(50)        :: wprob         ! archivo atomos
     integer(ip)          :: natoms        ! cantidad de atomos
     integer(ip)          :: nespec        ! cantidad especies
     integer(ip), pointer :: espe(:)       ! numero de especie
     integer(ip), pointer :: Tiespe(:)     ! tipo interno de especie
     integer(ip), pointer :: spin(:)       ! numero de spin por atomo
     integer(ip), pointer :: tipoPP(:)     ! tipo de PP
     real(rp),    pointer :: coorx(:)      ! pos x en order por atomo
     real(rp),    pointer :: coory(:)      ! pos y en order por atomo
     real(rp),    pointer :: coorz(:)      ! pos z en order por atomo
  end type atomo

  type especie
     character(50)        :: wprob1        ! name especie
     character(50)        :: wprob2        ! archivo ppseudo potencial
     integer(ip)          :: atnumber      ! numero atomico
     integer(ip)          :: ncp           ! nuemro cuantico ppal
     integer(ip)          :: ncl           ! nuemro cuantico l
     integer(ip)          :: ncm           ! nuemro cuantico m
     integer(ip)          :: nspin         ! nuemro spin
     integer(ip)          :: nlmax         ! orbital maximo
     integer(ip)          :: nlocal        ! orbital local
     real(rp)             :: valencia      ! electrones de valencia
     integer(ip)          :: nrad          ! numero de puntos que entrega el PP
     real(rp),allocatable :: radio(:)      ! coordenada radial del PP esferico
     real(rp),allocatable :: ppseu(:,:)    ! matriz dim nrad x nlmax con el PP
     real(rp),allocatable :: ppphi(:,:)    ! matriz dim nrad x nlmax con la PPPhi
     real(rp),allocatable :: nocupa(:)     ! ocupacion bien detalladita!
  end type especie

  !
  ! Chemical Species with their properties
  !
  type typ_speci
     character(5) :: name
     real(rp)     :: visco(2)
     integer(ip)  :: lawvi
     real(rp)     :: weigh
     real(rp)     :: densi(4)
     real(rp)     :: entha(2)
     real(rp)     :: prand
     real(rp)     :: lewis
     real(rp)     :: trang(5)
     real(rp)     :: cpcoe(10,4)
     real(rp)     :: activ(2)
  end type typ_speci

  !----------------------------------------------------------------------
  !
  ! Database for CFI combustion model
  !
  !----------------------------------------------------------------------

  type base_cfi

     integer(ip)   :: nclas           ! Number of particle classes
     integer(ip)   :: nccfi           ! Number of columns in table
     integer(ip)   :: nrcfi           ! Number of rows in table
     integer(ip)   :: nvcfi(5)        ! Number of values in table for each variable
     integer(ip)   :: nfcfi           ! Number of material properties, non-premixed
     integer(ip)   :: ndcfi           ! degree of freedom in tabulation
     real(rp)      :: imima(2)        ! minimum and maximum value of enthalpy
     real(rp)      :: fmima(2)        ! minimum and maximum value of mixture fraction mean
     real(rp)      :: inval(2,15)     ! inlet properties for non-premixed cases
     real(rp) , pointer :: ivcfi(:,:) ! subdivision of variables in table
     real(rp) , pointer :: table(:,:) ! Table for thermochemical database
     real(rp) , pointer :: ymass(:,:) ! Mass fraction at equilibrium for each mixture fraction

  end type base_cfi

  !----------------------------------------------------------------------
  !
  ! Mesh
  !
  !----------------------------------------------------------------------

  type cell

     real(rp)      :: coor(3,2)
     integer(ip)   :: neigh(6)
     integer(ip)   :: level
     integer(ip)   :: marked
     real(rp)      :: rsize

  end type cell

  !----------------------------------------------------------------------
  !
  ! Postprocess
  !
  !----------------------------------------------------------------------

  integer(ip), parameter ::        &
       nvars = 45,                 &    ! # set variables
       nvart = 10,                 &    ! # times for postprocess
       nvarw = 40,                 &    ! # witness point
       nvarp = 300                       ! # postprocess variables

  type typos
     integer(ip)           :: npp_inits                 ! Postprocess initial step
     integer(ip)           :: kfl_oonce(nvarp)          ! Postprocess only once
     integer(ip)           :: npp_stepi(nvarp)          ! Postprocess step interval for u,p, etc.
     integer(ip)           :: vox_stepi(nvarp)          ! Postprocess step interval for u,p, etc.
     integer(ip)           :: npp_iniso                 ! Postprocess initial condition
     integer(ip)           :: npp_stepw                 ! Postprocess witness interval
     integer(ip)           :: npp_setse(nvars)          ! Postprocess element sets calculation
     integer(ip)           :: npp_setsb(nvars)          ! Postprocess boundary sets calculation
     integer(ip)           :: npp_setsn(nvars)          ! Postprocess node sets calculation
     integer(ip)           :: per_setse(nvars)          ! Postprocess element sets calculation
     integer(ip)           :: per_setsb(nvars)          ! Postprocess boundary sets calculation
     integer(ip)           :: per_setsn(nvars)          ! Postprocess node sets calculation
     integer(ip)           :: npp_witne(nvarw)          ! Postprocess witness points
     integer(ip)           :: pos_alrea(nvarp)          ! Already postprocessed
     integer(ip)           :: vox_alrea(nvarp)          ! Already postprocessed
     integer(ip)           :: nfilt(nvarp)              ! Filter
     integer(ip)           :: nvaes                     ! Element set variables
     integer(ip)           :: nvabs                     ! Boundary set variables
     integer(ip)           :: nvans                     ! Node set variables
     integer(ip)           :: per_nvaes                 ! Element set variables to exchange
     integer(ip)           :: per_nvabs                 ! Boundary set variables to exchange
     integer(ip)           :: per_nvans                 ! Node set variables to exchange
     integer(ip)           :: nvawi                     ! Node witness variables
     integer(ip)           :: ipass                     ! Set memory allocated and header
     integer(ip)           :: lun_setse                 ! Element set unit imodu*10+6
     integer(ip)           :: lun_setsb                 ! Boundary set unit imodu*10+7
     integer(ip)           :: lun_setsn                 ! Node set unit imodu*10+8
     integer(ip)           :: lun_setsi                 ! Immersed bounday set
     integer(ip)           :: lun_witne                 ! Witness points
     character(150)        :: fil_setse
     character(150)        :: fil_setsb
     character(150)        :: fil_setsn
     character(150)        :: fil_setsi
     character(150)        :: fil_witne
     real(rp)              :: pos_tinit                 ! Postprocess initial time
     real(rp)              :: pos_times(nvart,nvarp)    ! Postprocess times for u,p, etc.
     real(rp)              :: pos_perio(nvarp)          ! Postprocess time period for u,p, etc.
     real(rp)              :: vox_times(nvart,nvarp)    ! Postprocess times for u,p, etc.
     real(rp)              :: vox_perio(nvarp)          ! Postprocess time period for u,p, etc.
     real(rp)              :: paese(5,nvars)            ! Element set parameters
     real(rp)              :: pabse(5,nvars)            ! Boundary set parameters
     real(rp)              :: panse(5,nvars)            ! Node set parameters
     character(5)          :: woese(nvars)              ! Name of the element set variables
     character(5)          :: wobse(nvars)              ! Name of the boundary set variables
     character(5)          :: wonse(nvars)              ! Name of the node set variables
     character(5)          :: wowit(nvarw)              ! Name of the witness
     character(5)          :: wopos(3,nvarp)            ! Name and character of postprocess variable
     real(rp),     pointer :: veset(:,:)                ! Set element values
     real(rp),     pointer :: vbset(:,:)                ! Set boundary values
     real(rp),     pointer :: vnset(:,:)                ! Set node values
     real(rp),     pointer :: viset(:,:)                ! Set IB values
     real(rp),     pointer :: witne(:,:)                ! Witness values
  end type typos

  !----------------------------------------------------------------------
  !
  ! Module
  !
  !----------------------------------------------------------------------

  type tymod
     type(typos),    pointer :: postp(:)           ! Postprocess
     type(soltyp),   pointer :: solve(:)           ! Algebraic solver
     type(soltyp),   pointer :: solad(:)           ! Adjoint algebraic solver
     type(eigtyp),   pointer :: eigen(:)           ! Eigenvalue solver
     type(bc_nodes), pointer :: tncod(:)           ! Node code type
     type(bc_nodes), pointer :: tgcod(:)           ! Geometrical node code type
     type(bc_bound), pointer :: tbcod(:)           ! Boundary code type
     integer(ip)             :: nvarn_bcs          ! Number of bc variables (nodes)
     integer(ip)             :: nvarg_bcs          ! Number of bc variables (geometrical)
     integer(ip)             :: nvarb_bcs          ! Number of bc variables (boundaries)
     integer(ip)             :: kfl_modul          ! Existence of module
     integer(ip)             :: kfl_delay          ! Delay module
     integer(ip)             :: kfl_conve          ! Cconvergence required
     integer(ip)             :: kfl_solve          ! When to solve module
     integer(ip)             :: kfl_goite          ! Keep on iterating
     integer(ip)             :: kfl_timei          ! Problem is transient
     integer(ip)             :: kfl_stead          ! Problem is steady
     integer(ip)             :: ndela              ! Steps to delay module
     integer(ip)             :: itinn              ! Module inner iteration
     integer(ip)             :: ittot              ! Total number of iteration
     integer(ip)             :: miinn              ! Total number of inner iteration
     integer(8)              :: mem_modul(2)       ! Module memory
     real(rp)                :: cpu_modul(30)      ! Module CPU time
     real(rp)                :: dtcri              ! Module critical time
     real(rp)                :: glres              ! Problem residuals
     character(6)            :: namod              ! Module name
     character(3)            :: exmod              ! Module extension
     integer(ip)             :: lun_pdata          ! File units
     integer(ip)             :: lun_outpu          ! ...
     integer(ip)             :: lun_conve
     integer(ip)             :: lun_rstpo
     integer(ip)             :: lun_rstar
     integer(ip)             :: lun_timin
     !integer(ip)             :: lun_time
     character(150)          :: fil_pdata          ! File names
     character(150)          :: fil_outpu          ! ...
     character(150)          :: fil_conve
     !character(150)          :: fil_time
     character(150)          :: fil_rstar
     character(150)          :: fil_timin
     type(comm_data_par), pointer :: commd         ! Communication type
  end type tymod

  !----------------------------------------------------------------------
  !
  ! Filters used in postprocess
  !
  !----------------------------------------------------------------------

  integer(ip), parameter   :: ncofi=3
  type tyfil
     integer(ip)           :: numbe              ! Filter number
     integer(ip)           :: ifilt(ncofi)       ! Filter number in module
     integer(ip)           :: modul(ncofi)       ! Module that computes it (0=master)
     real(rp)              :: param(ncofi,10)    ! Parameters
     integer(ip), pointer  :: lfilt(:)           ! List of filtered nodes
  end type tyfil

  !----------------------------------------------------------------------
  !
  ! Element bin
  !
  !----------------------------------------------------------------------

  type typ_element_bin
     real(rp)             :: comin(3)
     real(rp)             :: comax(3)
     integer(ip)          :: boxes(3)
     integer(ip), pointer :: bin_size(:,:,:)
     type(i1p),   pointer :: list_elements(:,:,:)
  end type typ_element_bin

  !----------------------------------------------------------------------
  !
  ! OMPSS_DOMAIN
  !
  !----------------------------------------------------------------------

  type ompss_domain
     integer(ip), pointer :: neighbours(:)
     integer(ip), pointer :: elements(:)
     integer(ip)          :: neighIdx
     integer(ip)          :: elemIdx
  end type ompss_domain

end module def_kintyp
!> @}
