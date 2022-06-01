module def_kintyp
  !-----------------------------------------------------------------------
  !****f* defmod/def_kintyp
  ! NAME
  !   def_kintyp
  ! DESCRIPTION
  !   Definition of kinds and types. 
  !   "The range of the default integers is not specified in the language
  !   but on a computer with a word size of n bits, is often from 
  !   -2^{n-1} to +2^{n-1}-1. Thus on a 32-bit computer the range is
  !   often -2.14*10^9 to +2.14*10^9."
  !   M. Metclaf and J. Reid, FORTRAN 90/95 explained, 2nd edition.
  !
  !   Defaults are:
  !   Integers: 4-bytes 
  !   Reals:    8-bytes
  !
  !***
  !-----------------------------------------------------------------------

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
  type i1pp
     integer(ip)          :: n
     integer(ip), pointer :: l(:)
  end type i1pp
  type i1pi1p
     type(i1p),   pointer :: l(:)
  end type i1pi1p

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
     integer(ip)              :: kfl_funno
     integer(ip)              :: kfl_fixrs
     character(5)             :: cotag
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
     character(5)             :: cotag
     real(rp),        pointer :: bvnat(:)
  end type bc_bound1
  type bc_bound 
     integer(ip)              :: ndofn     
     integer(ip)              :: ncode     
     type(bc_bound1), pointer :: l(:)
  end type bc_bound
 
  !----------------------------------------------------------------------
  !
  ! Element data base type
  !
  !----------------------------------------------------------------------

  type elm
     ! User integration rule
     real(rp),    pointer :: shape(:,:)
     real(rp),    pointer :: deriv(:,:,:)
     real(rp),    pointer :: heslo(:,:,:)
     real(rp),    pointer :: weigp(:)
     real(rp),    pointer :: shaga(:,:)
     real(rp),    pointer :: posgp(:,:)
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
  ! ADR eqn
  !
  !----------------------------------------------------------------------

  type adrtyp
     integer(ip) :: kfl_timei_adr               ! ADR eqn: Time flag  
     integer(ip) :: kfl_advec_adr               ! ADR eqn: Advection flag
     integer(ip) :: kfl_diffu_adr               ! ADR eqn: Diffusion flag
     integer(ip) :: kfl_react_adr               ! ADR eqn: Reaction flag
     integer(ip) :: kfl_grdif_adr               ! ADR eqn: Diffusion gradient flag
     integer(ip) :: kfl_tisch_adr               ! ADR eqn: Time scheme flag
     integer(ip) :: kfl_taust_adr               ! ADR eqn: Stab. strategy
     integer(ip) :: kfl_sgsti_adr               ! ADR eqn: SGS tracking in time
     integer(ip) :: kfl_sgsno_adr               ! ADR eqn: SGS tracking in non-linearity
     integer(ip) :: kfl_tiacc_adr               ! ADR eqn: Time scheme accuracy
     integer(ip) :: kfl_shock_adr               ! ADR eqn: shcok capturing flag
     integer(ip) :: kfl_ellen_adr               ! ADR eqn: element length strategy
     integer(ip) :: kfl_stead_adr               ! ADR eqn: steady state flag
     real(rp)    :: bemol_adr                   ! ADR eqn: Integration by parts convective term
     real(rp)    :: pabdf_adr(10)               ! ADR eqn: BDF coefficients
     real(rp)    :: staco_adr(4)                ! ADR eqn: stability constants
     real(rp)    :: shock_adr                   ! ADR eqn: shock capturing coefficients
     real(rp)    :: safet_adr                   ! ADR eqn: safety factor
  end type adrtyp

  !----------------------------------------------------------------------
  !
  ! Solver
  !
  !----------------------------------------------------------------------

  type soltyp
     character(50)        :: wprob              ! Names
     character(50)        :: wsolv 
     character(50)        :: wprec 
     integer(ip)          :: nprob
     integer(ip)          :: nzmat              ! System size
     integer(ip)          :: nzrhs
     integer(ip)          :: nequa                
     integer(ip)          :: nunkn                
     integer(ip)          :: ndofn
     integer(ip)          :: ndof2
     integer(ip)          :: nsist
     integer(ip)          :: nesky
     integer(ip)          :: nseqn
     integer(ip)          :: nrhss
     integer(ip)          :: kfl_algso          ! Solver parameters
     integer(ip)          :: kfl_symme
     integer(ip)          :: kfl_symeq
     integer(ip)          :: kfl_facto
     integer(ip)          :: kfl_resid          ! Residual for iterative solvers
     integer(ip)          :: kfl_recov          ! Recover original residual after slave exchange
     integer(ip)          :: kfl_ortho          ! Orthogonolization GMRES
     integer(ip)          :: kfl_limit          ! Algebraic limiter
     integer(ip)          :: kfl_posgr          ! Postprocess groups
     integer(ip)          :: kfl_adres          ! Automatic residual
     integer(ip)          :: kfl_cmplx          ! Real or complex solver JELENA
     integer(ip)          :: kfl_assem          ! Matrix has been assembled
     integer(ip)          :: kfl_schur          ! Schur solver
     integer(ip)          :: kfl_zones          ! Zone
     integer(ip)          :: kfl_blogs          ! Block Gauss-Seidel treatment
     integer(ip)          :: miter             
     integer(ip)          :: nkryd
     !
     ! Fixity
     !
     integer(ip)          :: kfl_iffix          ! If solver takes care of imposing Dirichlet b.c.
     real(rp),    pointer :: bvess(:,:)         ! Dirichlet value
     integer(ip), pointer :: kfl_fixno(:,:)     ! Fixity array
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
     integer(ip), pointer :: ia(:)              ! Deflated CG: Sparse graph
     integer(ip), pointer :: ja(:)              ! Deflated CG: Sparse graph
     integer(ip), pointer :: amgro(:,:,:)       ! Deflated CG: Sparse graph
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
     integer(ip), pointer :: iLdef(:)           ! Deflated CG: sparse coarse matrix
     integer(ip), pointer :: jLdef(:)           ! Deflated CG: sparse coarse matrix
     integer(ip), pointer :: iUdef(:)           ! Deflated CG: sparse coarse matrix
     integer(ip), pointer :: jUdef(:)           ! Deflated CG: sparse coarse matrix   
     real(rp),    pointer :: Lndef(:)           ! Deflated CG: sparse coarse matrix    
     real(rp),    pointer :: Undef(:)           ! Deflated CG: sparse coarse matrix     
     integer(ip), pointer :: invpRdef(:)        ! Deflated CG: sparse coarse matrix
     integer(ip), pointer :: invpCdef(:)        ! Deflated CG: sparse coarse matrix
     real(rp),    pointer :: askyldef(:)        ! Deflated CG: skyline coarse matrix   
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
     ! Penalization
     !
     integer(ip)          :: kfl_penal          ! Penalization: method
     real(rp)             :: penal              ! Penalization: parameter
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
     ! Sparse direct solver
     !
     integer(ip)          :: kfl_alcsr          ! If memory has been allocated
     integer(ip), pointer :: iL(:)
     integer(ip), pointer :: jL(:)
     integer(ip), pointer :: iU(:)
     integer(ip), pointer :: jU(:)    
     real(rp),    pointer :: Ln(:)    
     real(rp),    pointer :: Un(:)    
     !
     ! AII preconditioner
     !
     integer(ip), pointer :: iLpre(:)
     integer(ip), pointer :: jLpre(:)
     integer(ip), pointer :: iUpre(:)
     integer(ip), pointer :: jUpre(:)    
     real(rp),    pointer :: Lnpre(:)    
     real(rp),    pointer :: Unpre(:)      
     integer(ip), pointer :: invpR(:)
     integer(ip), pointer :: invpC(:)
     real(rp),    pointer :: Aiipre(:,:,:)         
     !
     ! Deflated multigrid preconditioner
     !
     integer(ip), pointer :: iLpredef(:)
     integer(ip), pointer :: jLpredef(:)
     integer(ip), pointer :: iUpredef(:)
     integer(ip), pointer :: jUpredef(:)    
     real(rp),    pointer :: Lnpredef(:)    
     real(rp),    pointer :: Unpredef(:)      
     integer(ip), pointer :: invpRpredef(:)
     integer(ip), pointer :: invpCpredef(:)
     real(rp),    pointer :: askylpredef(:)         
     integer(ip)          :: kfl_defpr          ! Smoother preconditioner
     !
     ! Coarse Aii
     !
     integer(ip), pointer :: iaaii(:)
     integer(ip), pointer :: jaaii(:)
     integer(ip), pointer :: iLaii(:)
     integer(ip), pointer :: jLaii(:)
     integer(ip), pointer :: iUaii(:)
     integer(ip), pointer :: jUaii(:)    
     real(rp),    pointer :: Lnaii(:)    
     real(rp),    pointer :: Unaii(:)      
     integer(ip), pointer :: invpRaii(:)
     integer(ip), pointer :: invpCaii(:)
     integer(ip), pointer :: iskylaii(:)        ! Coarse Aii: Skyline index
     integer(ip), pointer :: idiagaii(:)        ! Coarse Aii: Pointer to skyline diagonal
     real(rp),    pointer :: aiicoarse(:,:,:)         
     integer(ip), pointer :: lgaii(:)           ! List of group
     integer(ip)          :: kfl_deaii          ! Type of direct solver
     integer(ip)          :: ngaii              ! Number of groups
     !
     ! Other
     !
     real(rp)             :: solco
     real(rp)             :: relso
     real(rp)             :: dtinv   
     real(rp)             :: xdiag              ! Coefficient for diagonal solver
     real(rp)             :: resin              ! Initial preconditioned residual
     real(rp)             :: resfi              ! Final preconditioned residual
     real(rp)             :: resi2              ! Initial residual
     real(rp)             :: resf2              ! Final residual
     real(rp)             :: reni2              ! Non-normalized initial residual
     real(rp)             :: cputi(10)          ! CPU time
     real(rp)             :: cpu_schur(10)      ! Schur complement solver timing
     real(rp)             :: adres              ! Adaptive residual tolerance
     real(rp)             :: solmi              ! Minimum solver tolerance
     real(rp)             :: xorth              ! Orthogonality check
     integer(ip)          :: heade              ! Track if header has been written in file
     integer(ip)          :: itsol(3)           ! Solver statistics
     integer(ip)          :: nsolv
     integer(ip)          :: kfl_cvgso          ! Output
     integer(ip)          :: lun_cvgso
     character(150)       :: fil_cvgso
     integer(ip)          :: kfl_solve
     integer(ip)          :: lun_solve
     integer(ip)          :: kfl_exres          ! Exact residual output
     character(150)       :: fil_solve
     integer(ip)          :: kfl_preco          ! Preconditioner
     integer(ip)          :: kfl_leftr          ! Left(0) or right(1)
     integer(ip)          :: itpre              ! Preconditioner iterations
     integer(ip)          :: kfl_marke          ! Output of matrix in market format
     integer(ip)          :: kfl_force          ! Force continuity after calling solver
     integer(ip)          :: nzpre
     integer(ip)          :: lfill
     real(rp)             :: thres
     integer(ip), pointer :: lpdof(:)
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
     integer(ip)          :: kfl_facto
     integer(ip)          :: kfl_massm          ! Mass matrix type
     integer(ip)          :: miter              ! Solver number of iterations
     integer(ip)          :: itsol(3)           ! Solver statistics
     integer(ip)          :: nsolv
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
  !
  ! Parallelization 
  !
  type comm_data_par
     integer(ip)            :: nneig
     integer(ip), pointer   :: neights(:)
     ! Interior, own and other boundary nodes
     integer(ip)            :: npoi1
     integer(ip)            :: npoi2
     integer(ip)            :: npoi3
     integer(ip)            :: npoi4
     ! Node communication
     integer(ip), pointer   :: bound_size(:)
     integer(ip), pointer   :: bound_perm(:)
     integer(ip), pointer   :: bound_adja(:)
     integer(ip), pointer   :: bound_scal(:)
     integer(ip)            :: bound_dim
     ! Face communication
     integer(ip), pointer   :: bface_size(:)
     integer(ip), pointer   :: bface_perm(:)
     integer(ip)            :: bface_dim
     ! Fringe node communication
     integer(ip), pointer   :: frins_size(:)
     integer(ip), pointer   :: frins_perm(:)
     integer(ip)            :: frins_dim
     integer(ip), pointer   :: frinr_size(:)
     integer(ip), pointer   :: frinr_perm(:)
     integer(ip)            :: frinr_dim
     ! Communicator
     integer(4)             :: PAR_COMM_WORLD
  end type comm_data_par

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
       nvars = 25,                 &    ! # set variables
       nvart = 10,                 &    ! # times for postprocess
       nvarw = 20,                 &    ! # witness point
       nvarp = 70                       ! # postprocess variables

  type typos
     integer(ip)           :: npp_inits                 ! Postprocess initial step
     integer(ip)           :: kfl_oonce(nvarp)          ! Postprocess only once
     integer(ip)           :: npp_stepi(nvarp)          ! Postprocess step interval for u,p, etc.
     integer(ip)           :: vox_stepi(nvarp)          ! Postprocess step interval for u,p, etc.
     integer(ip)           :: npp_iniso                 ! Postprocess initial condition
     integer(ip)           :: npp_setse(nvars)          ! Postprocess element sets calculation
     integer(ip)           :: npp_setsb(nvars)          ! Postprocess boundary sets calculation
     integer(ip)           :: npp_setsn(nvars)          ! Postprocess node sets calculation
     integer(ip)           :: npp_witne(nvarw)          ! Postprocess witness points
     integer(ip)           :: pos_alrea(nvarp)          ! Already postprocessed
     integer(ip)           :: vox_alrea(nvarp)          ! Already postprocessed
     integer(ip)           :: nfilt(nvarp)              ! Filter
     integer(ip)           :: nvaes                     ! Element set variables
     integer(ip)           :: nvabs                     ! Boundary set variables
     integer(ip)           :: nvans                     ! Node set variables
     integer(ip)           :: nvawi                     ! Node witness variables
     integer(ip)           :: ipass                     ! Set memory allocated and header
     integer(ip)           :: lun_setse                 ! Element set unit imodu*10+6
     integer(ip)           :: lun_setsb                 ! Boundary set unit imodu*10+7
     integer(ip)           :: lun_setsn                 ! Node set unit imodu*10+8
     integer(ip)           :: lun_setsi    
     integer(ip)           :: lun_witne  
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
     type(eigtyp),   pointer :: eigen(:)           ! Eigenvalue solver
     type(bc_nodes), pointer :: tncod(:)           ! Node code type
     type(bc_nodes), pointer :: tgcod(:)           ! Geometrical node code type
     type(bc_bound), pointer :: tbcod(:)           ! Boundary code type
     integer(ip)             :: nvari_sol          ! Number of solver variables
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
     integer(ip)             :: miinn(3)           ! Total number of inner iteration
     integer(8)              :: mem_modul(2)       ! Module memory
     real(rp)                :: cpu_modul(30)      ! Module CPU time
     real(rp)                :: dtcri              ! Module critical time
     real(rp)                :: glres              ! Problem residuals
     character(6)            :: namod              ! Module name
     character(3)            :: exmod              ! Module extension
     integer(ip)             :: lun_pdata          ! File units
     integer(ip)             :: lun_outpu          ! ...
     integer(ip)             :: lun_conve 
     integer(ip)             :: lun_rstar 
     character(150)          :: fil_pdata          ! File names
     character(150)          :: fil_outpu          ! ...
     character(150)          :: fil_conve   
     character(150)          :: fil_rstar  
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
  ! Types of Lagrangian types of particles
  !
  !----------------------------------------------------------------------
  integer(ip), parameter  :: mlapr = 10        ! Maximum number of properties in Lagrangian particles
  type typlatyp
     integer(ip)            :: kfl_exist         ! If type exists
     integer(ip)            :: kfl_modla         ! Transport model
     integer(ip)            :: kfl_grafo         ! Gravity
     integer(ip)            :: kfl_buofo         ! Buoyancy
     integer(ip)            :: kfl_drafo         ! Drag force model
     integer(ip)            :: kfl_brown         ! Drag force model
     integer(ip)            :: kfl_extfo         ! External force
     real(rp)               :: denpa             ! Density particle
     real(rp)               :: spher             ! Particle sphericity
     real(rp)               :: diame             ! Particle diameter
     real(rp)               :: calor             ! Calorific capacity
     real(rp)               :: emisi             ! Particle radiation emissivity
     real(rp)               :: scatt             ! Particle radiation scattering factor
     real(rp)               :: diffu             ! Particle diffusion coefficient [m^2/s]
     real(rp)               :: prope(mlapr)      ! Default value of properties for each type
     integer(ip)            :: prova(mlapr)      ! ID of variable in properties vector
  end type typlatyp

end module def_kintyp
