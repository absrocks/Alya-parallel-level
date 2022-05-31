!===================================================================================!
!
!< 2014Sep24 -> commdom_plepp_on_field
!< 2014Nov04 -> commdom_plepp_set_vals
!< 2014Dic14 -> relax
!< 2015Jan29 -> commdom_plepp_set_source_nodes
!< 2015Feb18 -> commdom_plepp_check_fixno
!< 2015Feb18 -> commdom_plepp_set_vals + 'debug, norm2'
!< 2015Feb18 -> commdom_plepp_reduce_sum
!< 2015Abr02 ->
!< 2016Mar22 ->
!< 2016Mar23 ->
!< 2016MAR30 ->
!< 2017JAN07 -> (COMMDOM<=-3)||(COMMDOM>=3)
!
!===================================================================================!
module mod_commdom_plepp
  use mod_commdom_alya, only: INONE
#ifdef COMMDOM
  use def_parame, only: ip, rp
  use def_kermod, only: ndivi
  use def_master, only: inotmaster, imaster, isequen, islave, inotslave, iparall
  use def_master, only: routp, dtime, kfl_paral
  use def_master, only: ittim, cutim, mitim
  use def_master, only: tempe
  use def_domain, only: coord, mnode, nelem, ndime, npoin
  use def_domain, only: LESET !meshe, LBSET
  use def_domain, only: ltype, nnode, ngaus, lnods, coord
  use mod_parall, only: PAR_COMM_CURRENT, PAR_UNIVERSE_SIZE
  use mod_parall, only: PAR_COMM_WORLD, PAR_COMM_UNIVERSE, PAR_COMM_MY_CODE
  use mod_std
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip),   parameter :: cp = 64
  character(cp), parameter :: frmt = '(E11.4)'   !<-- frmt = '(f8.2)'
  integer(ip)              :: n_vertices_i, n_elements_i, n_vertices_j
  real(rp),    allocatable :: vec_ij(:), vec_ji(:,:)
  real(rp),    allocatable :: tetra_coords_j(:,:)
  real(rp),    allocatable :: prop_i(:), prop_j(:)

  type COMMDOM_PLEPP_COUPLING    !< type(COMMDOM_PLEPP_COUPLING) CPLNG
    integer(ip)              :: local_comm           = -1_ip
    integer(ip)              :: commij               = -1_ip
    integer(ip)              :: commji               = -1_ip
    integer(ip)              :: n_vertices_i         = -1_ip
    integer(ip)              :: n_elements_i         = -1_ip
    integer(ip)              :: n_vertices_j         = -1_ip
    integer(ip)              :: n_ij                 = -1_ip
    integer(ip)              :: n_ji                 = -1_ip
    integer(ip)              :: stride               = -1_ip
    !
    real(rp)             :: tol                  =  0.0_rp
    real(rp),    pointer :: var_ij(:,:)            => null()
    real(rp),    pointer :: var_ji(:,:)            => null()
    real(rp),    allocatable :: vec_ij(:,:)          !=> null()
    real(rp),    allocatable :: vec_ji(:,:)          !=> null()
    real(rp),    pointer :: tetra_coords_j(:,:)  => null()
    !
    real(rp),    pointer :: dist_coords_j(:)     => null()
    integer(ip), pointer :: dist_locations_i(:)  => null()
    integer(ip), pointer :: interior_list_j(:)   => null()
    real(rp),    pointer :: dist_props_j(:)      => null()

    integer(ip), pointer :: idx_coords_j(:)      => null()
    !
    character(cp)  :: namei       = ''
    character(cp)  :: namej       = ''
    character(cp)  :: app_type    = ''
    character(cp)  :: app_name    = ''
    character(cp)  :: module_name = ''
  end type COMMDOM_PLEPP_COUPLING
  !
  type(COMMDOM_PLEPP_COUPLING) PLEPP_CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  private
    interface commdom_plepp_exchange
      module procedure &
                       commdom_plepp_exchange01, &
                       commdom_plepp_exchange02
    end interface commdom_plepp_exchange
    public :: commdom_plepp_exchange

    public :: commdom_plepp_init
    public :: commdom_plepp_set_mesh
    public :: commdom_plepp_coupling_send_msg
    public :: commdom_plepp_exchange01
    public :: commdom_plepp_exchange02
    public :: commdom_plepp_compare_dtinv
    public :: commdom_plepp_set_vals
    public :: commdom_plepp_set_source_nodes
    public :: commdom_plepp_check_fixno
    public :: commdom_plepp_reduce_sum
    public :: PLEPP_CPLNG
    public :: commdom_plepp_on_field
    public :: commdom_plepp_locator_send_nodal_var00
    public :: commdom_plepp_locator_send_nodal_var01  !< 2015Abr02
    public :: COMMDOM_PLEPP_COUPLING
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  contains

!-------------------------------------------------------------------------||---!
!-----------------------------------------------------| par_split_universe |---!
!-------------------------------------------------------------------------||---!

  !============================================================================!
  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date
  !> @brief
  !> @details
  !>
  !>  @code
  !>  ./Sources/kernel/parall/par_code_split_universe.f90
  !>   ...
  !>   call par_commdom_init()      !> JMAKE
  !>   if( PAR_UNIVERSE_SIZE > 1 )
  !>   ...
  !>  @endcode
  !>
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_init()
  implicit none

  integer(ip)    :: world_comm
  character(cp)  :: send, recv = ''
  character(cp)  :: token
  integer(ip)    :: app_ok, n_send, n_recv
  integer(ip)    :: oki, okj
  integer(4)     :: iarg

  PLEPP_CPLNG%module_name =  ''

  !-----------------------------------------------------------------------||---!
  !app_type = "ALYA_CFD"  !< coupling with syrthes
  PLEPP_CPLNG%app_type = "SYRTHES 4" !< coupling with saturne

  PLEPP_CPLNG%namei = "FLUID"
  PLEPP_CPLNG%namej = "SOLID"
!#if COMMDOM>=3   !< 2017JAN07
#if (COMMDOM<=-3)||(COMMDOM>=3)
  PLEPP_CPLNG%namei = "DIRIC"        !< 2016Mar22
  PLEPP_CPLNG%namej = "NEUMA"        !< 2016Mar22
#endif

  PLEPP_CPLNG%tol   =  1.0e-3
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!

  call commdom_create()

  do iarg = 1, command_argument_count() !iargc()
    token = ""
    call get_command_argument(iarg, token)
    call commdom_set_argvs(trim(token), len_trim(token))
  enddo

  token = "--name"
  call commdom_analyse_argvs(trim(token), len_trim(token))

  PLEPP_CPLNG%app_name = ""
  call commdom_get_argvs(PLEPP_CPLNG%app_name)
  call commdom_set_names(trim(PLEPP_CPLNG%app_type), len_trim(PLEPP_CPLNG%app_type), trim(PLEPP_CPLNG%app_name), len_trim(PLEPP_CPLNG%app_name))
  PLEPP_CPLNG%module_name = trim(PLEPP_CPLNG%app_name)
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  world_comm = PAR_COMM_UNIVERSE
  call commdom_create_commij(world_comm, PLEPP_CPLNG%local_comm)

  PAR_COMM_WORLD   = PLEPP_CPLNG%local_comm
  PAR_COMM_CURRENT = PLEPP_CPLNG%local_comm
  !PAR_COMM_MY_CODE = local_comm
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  PLEPP_CPLNG%commij = -1
  PLEPP_CPLNG%commji = -1

  oki = -1
  okj = -1

  call commdom_who_iam(trim(PLEPP_CPLNG%namei), len_trim(PLEPP_CPLNG%namei), oki)
  if(oki==1) then
    call commdom_who_areyou(trim(PLEPP_CPLNG%namej), len_trim(PLEPP_CPLNG%namej), okj)
    if(okj==1) then
      call commdom_get_commij(trim(PLEPP_CPLNG%namej), len_trim(PLEPP_CPLNG%namej), PLEPP_CPLNG%commij)
    endif
  endif

  oki = -1
  okj = -1

  call commdom_who_iam(trim(PLEPP_CPLNG%namej), len_trim(PLEPP_CPLNG%namej), okj)
  if(okj==1) then
    call commdom_who_areyou(trim(PLEPP_CPLNG%namei), len_trim(PLEPP_CPLNG%namei), oki)
    if(oki==1) then
      call commdom_get_commij(trim(PLEPP_CPLNG%namei), len_trim(PLEPP_CPLNG%namei), PLEPP_CPLNG%commji)
    endif
  endif

  send = trim(PLEPP_CPLNG%app_name)
 !call commdom_plepp_coupling_send_msg( len_trim(send), trim(send), recv)      !< 2016MAR30
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
    PAR_UNIVERSE_SIZE = -1 !=> if( PAR_UNIVERSE_SIZE > 1 ) ...
  endif
  !-----------------------------------------------------------------------||---!

  end subroutine
  !=============================================================================!

!-------------------------------------------------------------------------||---!
!-----------------------------------------------------------------| domain |---!
!-------------------------------------------------------------------------||---!

  !=================================================================================!
  !---------------------------------------------------------------------------------!
  !> @author  JM Zavala-Ake
  !> @date
  !> @brief
  !> @details
  !---------------------------------------------------------------------------------!
  !> +
  !> |_./Sources/kernel/master/Turnon.f90
  !>   |_./Sources/kernel/domain/domain()
  !>     |_...
  !>     |_ if( IPARALL ) call renelm()
  !>     |_ call par_commdom_set_mesh02() !> JMAKE
  !>     |_ ...
  !>
  !---------------------------------------------------------------------------------!
  subroutine commdom_plepp_set_mesh( id_fixbo_j, stride )
  implicit none
  integer(ip), optional, intent(in) :: id_fixbo_j
  integer(ip), optional, intent(in) :: stride
  !
  integer(ip)   :: oki, okj
  integer(ip)   :: n_recv, n_send
  integer(ip)   :: n_dime_i, n_node_i
  !
  real(rp),    pointer :: vertex_coords_i(:,:) => null()
  integer(ip), pointer ::    vertex_num_i(:,:) => null()
  integer(ip), pointer ::   vertex_type_i(:  ) => null()
  real(rp),    pointer :: vertex_coords_j(:,:) => null()
  integer(ip), pointer :: interior_list_j(:  ) => null()
  !
  PLEPP_CPLNG%n_ij   = -1
  PLEPP_CPLNG%n_ji   = -1
  PLEPP_CPLNG%stride =  1
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    n_vertices_i = 0
    n_elements_i = 0
    n_vertices_j = 0
    n_dime_i = 0
    n_node_i = 0
    if(INOTMASTER) then
      n_dime_i = ndime
      n_node_i = mnode
      !
      n_vertices_i    =  npoin
      n_elements_i    =  nelem
      vertex_coords_i => coord
      vertex_num_i    => lnods
      vertex_type_i   => ltype

      n_vertices_j    =  npoin
      vertex_coords_j => coord
    else
      allocate( vertex_coords_i(0,0) )
      allocate(    vertex_num_i(0,0) )
      allocate(   vertex_type_i(0  ) )
      allocate( vertex_coords_j(0,0) )
    endif
    !
    if( present(id_fixbo_j) ) then
      if(id_fixbo_j>0) call commdom_plepp_on_field(PLEPP_CPLNG, id_fixbo_j, n_vertices_j, vertex_coords_j)
    endif
    !---------------------------------------------------------------------||---!

    !---------------------------------------------------------------------||---!
    if( .not.(PLEPP_CPLNG%commij == -1) ) then
      call commdom_locator_create2(PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commij, PLEPP_CPLNG%tol)
    endif
    if( .not.(PLEPP_CPLNG%commji == -1) ) then
      call commdom_locator_create2(PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commji, PLEPP_CPLNG%tol)
    endif

!
!    call commdom_locator_set_mesh(   n_vertices_i, &
!                                     n_elements_i, &
!                                  vertex_coords_i, &
!                                     vertex_num_i, &
!                                     n_vertices_j, &
!                                  vertex_coords_j )
!
    call commdom_locator_set_cs_mesh(   n_vertices_i, &
                                        n_elements_i, &
                                     vertex_coords_i, &
                                        vertex_num_i, &
                                       vertex_type_i, &
                                        n_vertices_j, &
                                     vertex_coords_j, ndime )  !< 2016Mar22

    !< intel -g fail
!    call commdom_locator_set_mesh(                     n_vertices_i, &
!                                                       n_elements_i, &
!                                    meshe(0)%coord(1:ndime,1:npoin), &
!                                    meshe(0)%lnods(1:mnode,1:nelem), &
!                                                       n_vertices_j, &
!                                    meshe(0)%coord(1:ndime,1:npoin) )
    !
    call commdom_locator_save_dist_coords(0, PLEPP_CPLNG%local_comm)
    !
    !call commdom_plepp_coupling_edf_start()
    !---------------------------------------------------------------------||---!

    !---------------------------------------------------------------------||---!
    n_recv = 0
    n_send = 0
    call commdom_locator_get_n_dist_points( n_send )
    call commdom_locator_get_n_interior(    n_recv )
    !
    PLEPP_CPLNG%n_ij   = n_send
    PLEPP_CPLNG%n_ji   = n_recv
    if( present(stride) ) PLEPP_CPLNG%stride = stride
    !
    if( .not.allocated(  tetra_coords_j)               ) allocate(               tetra_coords_j(          n_send, mnode) )
    if( .not.associated( PLEPP_CPLNG%dist_locations_i) ) allocate( PLEPP_CPLNG%dist_locations_i(          n_send       ) )
    if( .not.associated( PLEPP_CPLNG%dist_coords_j)    ) allocate( PLEPP_CPLNG%dist_coords_j(             n_send*3     ) )
    if( .not.associated( PLEPP_CPLNG%var_ij)           ) allocate( PLEPP_CPLNG%var_ij(PLEPP_CPLNG%stride, n_send       ) )
    if( .not.associated( PLEPP_CPLNG%var_ji)           ) allocate( PLEPP_CPLNG%var_ji(PLEPP_CPLNG%stride, n_recv       ) )
    if( .not.associated( PLEPP_CPLNG%interior_list_j)  ) allocate( PLEPP_CPLNG%interior_list_j(           n_recv       ) )
    !
    if(INOTMASTER) then
      PLEPP_CPLNG%dist_locations_i = -1
      PLEPP_CPLNG%dist_coords_j    = -1
      call commdom_locator_get_dist_locations( PLEPP_CPLNG%dist_locations_i )
      call commdom_locator_get_dist_coords(    PLEPP_CPLNG%dist_coords_j    )
      call commdom_plepp_locator_send_nodal_var00(PLEPP_CPLNG%dist_locations_i, PLEPP_CPLNG%dist_coords_j, n_send, tetra_coords_j)

!      PLEPP_CPLNG%interior_list_j(1:n_recv) = -1
!      call commdom_locator_get_interior_list( PLEPP_CPLNG%interior_list_j(1:n_recv) )

      if( present(id_fixbo_j) ) then
        if(id_fixbo_j>0) then
          allocate( interior_list_j( n_recv  ) )
          interior_list_j(1:n_recv) = -1
          call commdom_locator_get_interior_list( interior_list_j(1:n_recv) )
          PLEPP_CPLNG%interior_list_j(1:n_recv) = PLEPP_CPLNG%idx_coords_j( interior_list_j(1:n_recv)  )
          deallocate( interior_list_j )
        endif
      else
          PLEPP_CPLNG%interior_list_j(1:n_recv) = -1
          call commdom_locator_get_interior_list( PLEPP_CPLNG%interior_list_j(1:n_recv) )
      endif
!
!      leset = 0
!      LESET = 0
!      LESET(PLEPP_CPLNG%dist_locations_i) = -1
!      !print *, "min, nelem, max", minval(dist_locations_i), nelem, maxval(dist_locations_i), shape( LESET )
!
    endif
    !
  endif
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!

  end subroutine
  !=================================================================================!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  ! FIELDS, NUMBE=1
  !   FIELD=1, DIMENSION=1, NODE
  !     INCLUDE  ../Mesh01_CHARACTERISTICS.alya
  !   END_FIELD
  ! END_FIELDS
  !
  ! kfl_field(:,1:nfiel) ! DIM, NBOUN_TYPE, EXCLU
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_on_field(CPLNG, where_number, n_coords_j, coords_j)
  use def_domain, only:  nfiel, kfl_field, xfiel
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  type(COMMDOM_PLEPP_COUPLING), intent(inout)   :: CPLNG
  integer(ip),                  intent(in   )   :: where_number
  integer(ip),                  intent(inout  ) :: n_coords_j
  real(rp),    pointer,         intent(inout  ) :: coords_j(:,:)

  integer(ip) :: ipoin, ifiel, idime, i_coords_j
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(nfiel==0) then
    print *, "[commdom_plepp_on_field]", " Into '*.dom.dat' include -> 'FIELD=1, DIMENSION=1, NODE'"
    call runend('EXIT!!')
  endif
  !
  ifiel = 1_ip
  idime = 1_ip
  CPLNG%n_vertices_j = 0_ip
  if(INOTMASTER) then
    !
    do ipoin =1_ip,npoin
      if( int(xfiel(ifiel) % a(idime,ipoin,1)) == where_number) then
!        print*, coord(:,ipoin), xfiel(ifiel) % a(idime,ipoin)
        CPLNG%n_vertices_j = CPLNG%n_vertices_j + 1_ip
      endif
    enddo
    !
    !if(CPLNG%n_vertices_j==0) then
    !    print *, "n_where_number:"
    !    call runend('EXIT!!')
    !endif
    !
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(.not.associated(CPLNG%idx_coords_j)) allocate( CPLNG%idx_coords_j(CPLNG%n_vertices_j) )
  !
  if(INOTMASTER) then
    if(associated(coords_j)) then
      nullify(coords_j)
      allocate(coords_j(ndime,CPLNG%n_vertices_j) )
    endif
  endif
  !
  n_coords_j = CPLNG%n_vertices_j
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER.or.ISEQUEN) then
    i_coords_j = 0_ip
    do ipoin =1_ip,npoin
      if( int(xfiel(ifiel) % a(idime,ipoin,1)) == where_number) then
        i_coords_j = i_coords_j + 1_ip
                  coords_j(1:ndime,i_coords_j) = coord(1:ndime,ipoin)
        CPLNG%idx_coords_j(        i_coords_j) = ipoin
     endif
    enddo
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER.or.ISEQUEN) print *, "[commdom_plepp_on_field]", n_coords_j
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !=================================================================================!
  subroutine commdom_plepp_memall_init(CPLNG)
  use def_domain, only: nelem, ndime, npoin, nnode, ngaus
  use def_master, only: mem_modul
  use def_master, only: modul
  use mod_memory, only: memory_alloca
  implicit none
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( (CPLNG%commij /= -1).or.(CPLNG%commji /= -1) ) then
    call memory_alloca( mem_modul(1:2,modul),   'tetra_coords_j', 'commdom_plepp_memall',   CPLNG%tetra_coords_j, CPLNG%n_ij, mnode)
    call memory_alloca( mem_modul(1:2,modul), 'dist_locations_i', 'commdom_plepp_memall', CPLNG%dist_locations_i, CPLNG%n_ij)
    call memory_alloca( mem_modul(1:2,modul),    'dist_coords_j', 'commdom_plepp_memall',    CPLNG%dist_coords_j, CPLNG%n_ij*3_ip)
    call memory_alloca( mem_modul(1:2,modul),           'var_ij', 'commdom_plepp_memall',           CPLNG%var_ij, PLEPP_CPLNG%stride, CPLNG%n_ij)
    call memory_alloca( mem_modul(1:2,modul),           'var_ji', 'commdom_plepp_memall',           CPLNG%var_ji, PLEPP_CPLNG%stride, CPLNG%n_ji)
    call memory_alloca( mem_modul(1:2,modul),  'interior_list_j', 'commdom_plepp_memall',  CPLNG%interior_list_j, CPLNG%n_ji)
!
!    if( .not.allocated(  tetra_coords_j) ) allocate(  tetra_coords_j(  n_send,4) )
!    if( .not.allocated(dist_locations_i) ) allocate( dist_locations_i( n_send  ) )
!    if( .not.allocated(   dist_coords_j) ) allocate(    dist_coords_j( n_send*3) )
!    if( .not.allocated(          var_ij) ) allocate(           var_ij( n_send  ) )
!    if( .not.allocated(          var_ji) ) allocate(           var_ji( n_recv  ) )
!    if( .not.allocated( interior_list_j) ) allocate(  interior_list_j( n_recv  ) )
!
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !=================================================================================!

  !=================================================================================!
  subroutine commdom_plepp_memall_end(CPLNG)
  implicit none
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !=================================================================================!

!-------------------------------------------------------------------------||---!
!------------------------------------------------------------------| texts |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_check_fixno(fixno, idofn, fixval, ToDo)
  use def_domain, only: kfl_codno
  implicit none
  integer(ip),  intent(inout), pointer  :: fixno(:,:)
  integer(ip),  intent(in   )           :: idofn
  integer(ip),  intent(in   )           :: fixval
  logical(ip),  intent(in   )           :: ToDo
  !
  integer(ip) :: n_fixval, n_coupling
  integer(ip) :: icplng, ipoin, ii
  !
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
    n_fixval   = 0
    n_coupling = PLEPP_CPLNG%n_ji
    if(INOTMASTER.or.ISEQUEN) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    !
    ! n_fixval: cuantos 'fixno' iguales a 'fixval' existen?
    !
    n_fixval = count( fixno(idofn,PLEPP_CPLNG%interior_list_j) == fixval ,KIND=ip)
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(n_fixval==n_coupling) then
      !-------------------------------------------------------------------||---!
      !-------------------------------------------------------------------||---!
    else
      !-------------------------------------------------------------------||---!
      if(ToDo) then
        ii = 0
        do icplng = 1,n_coupling
          ipoin = PLEPP_CPLNG%interior_list_j(icplng)
          if( fixno(idofn,ipoin) /= fixval ) then
            ii = ii+ 1
            print *, ii, ",", ipoin, ",", fixno(idofn,ipoin), ",", kfl_codno(:,ipoin)
          endif
        enddo
        if(PLEPP_CPLNG%commij /= -1) then
          print *, "[commdom_plepp_check_fixno_i] ", "[i, ipoin, FIXNO, codno]"
          print *, "[commdom_plepp_check_fixno_i] ","n_coupling-n_fixno:", PLEPP_CPLNG%n_ji-n_fixval, ", USE FIXNO->", fixval
        else if(PLEPP_CPLNG%commji /= -1) then
          print *, "[commdom_plepp_check_fixno_j] ", "[j, ipoin, FIXNO, codno]"
          print *, "[commdom_plepp_check_fixno_j] ","n_coupling-n_fixno:", PLEPP_CPLNG%n_ji-n_fixval, ", USE FIXNO->", fixval
        endif
        call runend('EXIT!!')
      endif
      !-------------------------------------------------------------------||---!
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    endif
  endif
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_set_source_nodes(prop_out)
  implicit none
  logical(ip), intent(out)  :: prop_out(npoin)
  !
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER.or.ISEQUEN) then
      prop_out(1_ip:npoin) = .false.
      prop_out(PLEPP_CPLNG%interior_list_j) = .true.
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_set_vals(prop_in, prop_out, relax_in, debug, norm2)
  implicit none
  real(rp), intent( in)  :: prop_in(npoin)
  real(rp), intent(out)  :: prop_out(npoin)
  !
  logical(ip), optional, intent(in   ) :: debug
  real(rp),    optional, intent(in   ) :: relax_in
  real(rp),    optional, intent(inout) :: norm2
  !
  real(rp) :: relax = 1.0
  if( present(relax_in) ) relax = relax_in
  !
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER.or.ISEQUEN) then
     !prop_out(interior_list_j) = prop_in(interior_list_j)
     !prop_out(PLEPP_CPLNG%interior_list_j) = prop_in(PLEPP_CPLNG%interior_list_j)
      !
      !< 2014Dic14
      ! u(n+1) = relax * u(n+1)  + (1.0-relax) * u(n)
      !    OUT = relax * OUT     + (1.0-relax) * IN
      prop_out(PLEPP_CPLNG%interior_list_j) =  prop_in(PLEPP_CPLNG%interior_list_j) * (    relax) + &
                                              prop_out(PLEPP_CPLNG%interior_list_j) * (1.0-relax)
      !-------------------------------------------------------------------||---!
      !                                                                        !
      !-------------------------------------------------------------------||---!
      if( present(debug).and.debug ) then
        print *, "[commdom_plepp_set_vals]", minval(  prop_in(PLEPP_CPLNG%interior_list_j) ), &
                                                sum(  prop_in(PLEPP_CPLNG%interior_list_j) )/PLEPP_CPLNG%n_ji, &
                                                maxval( prop_in(PLEPP_CPLNG%interior_list_j) ), "<--"
        print *, "[commdom_plepp_set_vals]", minval( prop_out(PLEPP_CPLNG%interior_list_j) ), &
                                                sum( prop_out(PLEPP_CPLNG%interior_list_j) )/PLEPP_CPLNG%n_ji, &
                                             maxval( prop_out(PLEPP_CPLNG%interior_list_j) ), "-->"
      endif
      !-------------------------------------------------------------------||---!
      !                                                                        !
      !-------------------------------------------------------------------||---!
      if( present(norm2) ) then
        norm2 = dot_product( prop_in(PLEPP_CPLNG%interior_list_j) - prop_out(PLEPP_CPLNG%interior_list_j), &
                             prop_in(PLEPP_CPLNG%interior_list_j) - prop_out(PLEPP_CPLNG%interior_list_j)  )
      endif
      !-------------------------------------------------------------------||---!
      !                                                                        !
      !-------------------------------------------------------------------||---!
    endif
    !---------------------------------------------------------------------||---!
  endif
  !
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
!  subroutine commdom_plepp_set_vals(prop_in, prop_out, relax_in)
!  implicit none
!  real(rp), intent( in)  :: prop_in(npoin)
!  real(rp), intent(out)  :: prop_out(npoin)
!  !
!  real(rp), optional, intent(in) :: relax_in
!  real(rp) :: relax = 1.0
!  if( present(relax_in) ) relax = relax_in
!  !
!  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
!    !---------------------------------------------------------------------||---!
!    !                                                                          !
!    !---------------------------------------------------------------------||---!
!    if(INOTMASTER.or.ISEQUEN) then
!     !prop_out(interior_list_j) = prop_in(interior_list_j)
!     !prop_out(PLEPP_CPLNG%interior_list_j) = prop_in(PLEPP_CPLNG%interior_list_j)
!!
!      !
!      !< 2014Dic14
!      ! u(n+1) = relax * u(n+1)  + (1.0-relax) * u(n)
!      !    OUT = relax * OUT     + (1.0-relax) * IN
!      prop_out(PLEPP_CPLNG%interior_list_j) =  prop_in(PLEPP_CPLNG%interior_list_j) * (    relax) + &
!                                              prop_out(PLEPP_CPLNG%interior_list_j) * (1.0-relax)
!!
!    endif
!    !---------------------------------------------------------------------||---!
!  endif
!  !
!  end subroutine
!  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_exchange01(prop_i, prop_j)
  implicit none
  real(rp),    intent( in)  :: prop_i(npoin)
  real(rp),    intent(out)  :: prop_j(npoin)
  real(rp)      :: daux(3)
  character(cp) :: saux(3)
  !
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
print *, "USE 'commdom_plepp_exchange01'"
call runend('EXIT!!')
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    PLEPP_CPLNG%var_ij(1_ip,1_ip:PLEPP_CPLNG%n_ij) = -1
    call commdom_plepp_locator_send_nodal_var01(prop_i(1_ip:npoin),                                                 & !<---
                                                PLEPP_CPLNG%dist_locations_i(1:PLEPP_CPLNG%n_ij), PLEPP_CPLNG%n_ij, &
                                                tetra_coords_j(1_ip:PLEPP_CPLNG%n_ij,1_ip:mnode),                   & !--->
                                                PLEPP_CPLNG%var_ij(1_ip,1_ip:PLEPP_CPLNG%n_ij))
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
!    if((INOTMASTER.or.ISEQUEN).and.(.not.PLEPP_CPLNG%commij==-1)) then
    if(INOTMASTER.or.ISEQUEN) then
      daux(1) = minval( PLEPP_CPLNG%var_ij(1_ip,1_ip:PLEPP_CPLNG%n_ij) )
      daux(2) =    sum( PLEPP_CPLNG%var_ij(1_ip,1_ip:PLEPP_CPLNG%n_ij) )/PLEPP_CPLNG%n_ij
      daux(3) = maxval( PLEPP_CPLNG%var_ij(1_ip,1_ip:PLEPP_CPLNG%n_ij) )
      !
      write(saux(1), frmt) daux(1)
      write(saux(2), frmt) daux(2)
      write(saux(3), frmt) daux(3)
      !
      !print *, "prop_i", minval(prop_i), sum(prop_i)/npoin, maxval(prop_i), PLEPP_CPLNG%n_ij
      if(PLEPP_CPLNG%n_ij/=0) print *, "["//trim(PLEPP_CPLNG%module_name)//"]->["//trim(saux(1))//trim(saux(2))//trim(saux(3))//"]"
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!

    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    PLEPP_CPLNG%var_ji(1_ip,1_ip:PLEPP_CPLNG%n_ji) = -1
    call commdom_locator_exchange_double_scalar(PLEPP_CPLNG%var_ij, PLEPP_CPLNG%var_ji) !> send, recv
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
!    if((INOTMASTER.or.ISEQUEN).and.(.not.PLEPP_CPLNG%commji==-1)) then
    if(INOTMASTER.or.ISEQUEN) then
      prop_j(PLEPP_CPLNG%interior_list_j) = PLEPP_CPLNG%var_ji(1_ip,1:PLEPP_CPLNG%n_ji)
      !
      daux(1) = minval( prop_j(PLEPP_CPLNG%interior_list_j) )
      daux(2) =    sum( prop_j(PLEPP_CPLNG%interior_list_j) )/PLEPP_CPLNG%n_ji
      daux(3) = maxval( prop_j(PLEPP_CPLNG%interior_list_j) )
      !
      write(saux(1), frmt) daux(1)
      write(saux(2), frmt) daux(2)
      write(saux(3), frmt) daux(3)
      !
      if(PLEPP_CPLNG%n_ji/=0) print *, "["//trim(PLEPP_CPLNG%module_name)//"]<-["//trim(saux(1))//trim(saux(2))//trim(saux(3))//"]"
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_exchange02(prop_i, prop_j, stride)
  implicit none
  real(rp),    intent( in) :: prop_i(:,:)
  real(rp),    intent(out) :: prop_j(:,:)
  integer(ip), intent( in) :: stride
  real(rp)      :: daux(4)
  character(cp) :: saux(4)
  integer(ip)   :: ii
  !
  !if(present(stride_in)) stride = stride_in
  !
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if((stride<1_ip).or.(stride>PLEPP_CPLNG%stride)) then
      print *, "[commdom_plepp_exchange02] ERROR: stride<PLEPP_CPLNG%stride:", stride, "<",PLEPP_CPLNG%stride
      call runend('EXIT!!')
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    PLEPP_CPLNG%var_ij(1_ip:stride,1_ip:PLEPP_CPLNG%n_ij) = -1

  if(INOTMASTER.or.ISEQUEN) then
    do ii = 1,stride
      call commdom_plepp_locator_send_nodal_var01(prop_i(ii,1:npoin),                                                 & !<---
                                                  PLEPP_CPLNG%dist_locations_i(1:PLEPP_CPLNG%n_ij), PLEPP_CPLNG%n_ij, &
                                                  tetra_coords_j(1_ip:PLEPP_CPLNG%n_ij,1_ip:mnode),                   & !--->
                                                  PLEPP_CPLNG%var_ij(ii,1_ip:PLEPP_CPLNG%n_ij))
      !---------------------------------------------------------------------||---!
      !                                                                          !
      !---------------------------------------------------------------------||---!
        daux(1) = minval( PLEPP_CPLNG%var_ij(ii,1_ip:PLEPP_CPLNG%n_ij) )
        daux(2) =    sum( PLEPP_CPLNG%var_ij(ii,1_ip:PLEPP_CPLNG%n_ij) )/PLEPP_CPLNG%n_ij
        daux(3) = maxval( PLEPP_CPLNG%var_ij(ii,1_ip:PLEPP_CPLNG%n_ij) )
        !
        write(saux(1), frmt) daux(1)
        write(saux(2), frmt) daux(2)
        write(saux(3), frmt) daux(3)
        write(saux(4), '(I2)') ii
        !
        if(PLEPP_CPLNG%n_ij/=0)print*,"["//trim(PLEPP_CPLNG%module_name)//trim(saux(4))//"]->["//trim(saux(1))//","//trim(saux(2))//","//trim(saux(3))//"]"
    enddo
  endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!

    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    PLEPP_CPLNG%var_ji(1_ip:stride,1_ip:PLEPP_CPLNG%n_ji) = -1
    call commdom_locator_exchange_double_stride(PLEPP_CPLNG%var_ij, PLEPP_CPLNG%var_ji, stride) !> send, recv
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!

    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER.or.ISEQUEN) then
    do ii = 1,stride
      prop_j(ii,PLEPP_CPLNG%interior_list_j) = PLEPP_CPLNG%var_ji(ii,1:PLEPP_CPLNG%n_ji)
      !
      daux(1) = minval( prop_j(ii,PLEPP_CPLNG%interior_list_j) )
      daux(2) =    sum( prop_j(ii,PLEPP_CPLNG%interior_list_j) )/PLEPP_CPLNG%n_ji
      daux(3) = maxval( prop_j(ii,PLEPP_CPLNG%interior_list_j) )
      !
      write(saux(1), frmt) daux(1)
      write(saux(2), frmt) daux(2)
      write(saux(3), frmt) daux(3)
      write(saux(4), '(I2)') ii
      !
      if(PLEPP_CPLNG%n_ji/=0)print *, "["//trim(PLEPP_CPLNG%module_name)//trim(saux(4))//"]<-["//trim(saux(1))//","//trim(saux(2))//","//trim(saux(3))//"]"
    enddo
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_reduce_sum(prop_in, prop_out)
  implicit none
  real(rp), intent( in)  :: prop_in
  real(rp), intent(out)  :: prop_out
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER.or.ISEQUEN) then
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(.not.(PLEPP_CPLNG%commij==-1)) then
      call commdom_reduce_sum_real(prop_in, prop_out, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commij)
    endif
    if(.not.(PLEPP_CPLNG%commji==-1)) then
      call commdom_reduce_sum_real(prop_in, prop_out, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commji)
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !=================================================================================!

!-------------------------------------------------------------------------||---!
!------------------------------------------------------------------| texts |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_coupling_send_msg(n_char, char_send, char_recv)
  implicit none
  character(*), intent(in)     :: char_send
  character(*), intent(inout)  :: char_recv
  integer(ip),  intent(in)     :: n_char
  character(cp) :: cdummy
  integer(ip)   :: cdummy_len
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
    cdummy     = ''
    cdummy_len = 0
    if(.not.(PLEPP_CPLNG%commij==-1)) then
      call commdom_sendrecv_char(char_send,     n_char,    cdummy, cdummy_len, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commij)
      if(IMASTER.or.ISEQUEN) print *, "['", trim(PLEPP_CPLNG%namei), "'->'",  trim(char_send), "']"
    endif
    if(.not.(PLEPP_CPLNG%commji==-1)) then
      call commdom_sendrecv_char(   cdummy, cdummy_len, char_recv,     n_char, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commji)
      if(IMASTER.or.ISEQUEN) print *, "['", trim(PLEPP_CPLNG%namej), "'<-'",  trim(char_recv), "']"
    endif
    !if(IMASTER.or.ISEQUEN) then
    !  print *, "[", trim(char_send), "<-", trim(char_recv), "]"
    !endif
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !=================================================================================!

  subroutine commdom_plepp_compare_dtinv(dt_inv)
    implicit none
    real(rp), intent(inout) :: dt_inv
    real(rp)      :: send, recv, dtmin
    integer(ip)   :: comm, n_send, n_recv
    character(64) :: char_send, char_recv
    character(16) :: str(3)
  !-----------------------------------------------------------------------||---!
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then

    !---------------------------------------------------------------------||---!
    if( .not.(PLEPP_CPLNG%commij == -1) ) then
      comm = PLEPP_CPLNG%commij
    endif
    if( .not.(PLEPP_CPLNG%commji == -1) ) then
      comm = PLEPP_CPLNG%commji
    endif
    !---------------------------------------------------------------------||---!

    send   =  dt_inv
    recv   = -1.0_rp
    n_send =  1
    n_recv =  1
    call commdom_sendrecv_real(send, n_send, recv, n_recv, PLEPP_CPLNG%local_comm, comm)
    call commdom_bcast_real(recv, n_recv, PLEPP_CPLNG%local_comm, comm)

    dtmin  = min(1.0_rp/send, 1.0_rp/recv)
    dt_inv = 1.0_rp/dtmin

    !---------------------------------------------------------------------||---!
    if(IMASTER.or.ISEQUEN) then
      write(str(1),'(e11.4)') 1.0_rp/send
      write(str(2),'(e11.4)') 1.0_rp/recv
      write(str(3),'(e11.4)') dtmin
      print *, ""
      print *, '[commdom_plepp_compare_dtinv] '// trim(PLEPP_CPLNG%module_name) //' '&
                //'['//trim(str(1))&
                //','//trim(str(2))&
                //','//trim(str(3))&
                //']'
    endif
    !---------------------------------------------------------------------||---!

  endif
  !-----------------------------------------------------------------------||---!
  end subroutine

  !---------------------------------------------------------------------------<!
  !------------------------------------------------------| INTERPOLATIONs |---<!
  !---------------------------------------------------------------------------<!

  !=================================================================================!
  !---------------------------------------------------------------------------------!
  !> @author  JM Zavala-Ake
  !> @date
  !> @brief
  !> @details
  !---------------------------------------------------------------------------------!
  subroutine commdom_plepp_locator_send_nodal_var00(elemts_i, coords_j, n_dist_j, tetracoords_j)
  implicit none
  integer(ip), intent( in) :: n_dist_j
  integer(ip), intent( in) :: elemts_i(n_dist_j)
  real(rp),    intent( in) :: coords_j(n_dist_j*ndime) !< coords_j(n_dist_j*3)
  real(rp),    intent(out) :: tetracoords_j(n_dist_j,mnode)

  integer(ip) :: ii, ielem, pelty, pnode
  integer(ip) :: vertices_i(mnode)
  real(rp)    :: vol_coords_j(mnode)
  real(rp)    :: point_j(ndime) !< point_j(3)
!  real(rp)    :: prop_j(4)
!  integer(ip) :: alya_type
!  real(rp)    :: tol
  !-----------------------------------------------------------------------||---!
  if(IMASTER) then
  else
    !---------------------------------------------------------------------||---!
    elementary_loop: &
    do ii = 1,n_dist_j
      ielem = elemts_i(ii)
      pelty = ltype(ielem)
      !-------------------------------------------------------------------||---!
      pnode = nnode(pelty)
      vertices_i(1:mnode) = -1
      vertices_i(1:pnode) = lnods(1:pnode,ielem)
!
!<    point_j(1:3) = coords_j(ii*3-2:ii*3)  !> (ii-1)*3+jj = ii*3 + [-2,-1,0]
!
!<    call commdom_locator_interpolation(pelty, PLEPP_CPLNG%tol, coord(1:ndime,1:npoin), vertices_i(1:pnode), point_j(1:3), vol_coords_j(1:pnode) )
!<   !call commdom_locator_tetra_interpolation( coord(1:ndime,1:npoin), vertices_i(1:4), point_j(1:3), vol_coords_j(1:4) )
!
      point_j(1:ndime) = coords_j(ndime*ii-ndime+1:ndime*ii)  !> ndime(ii-1)+[1:ndime] = [ ndime(ii-1)+1:ndime*ii ]                                      !< 2016Mar23
      call commdom_locator_interpolation(pelty, PLEPP_CPLNG%tol, coord(1:ndime,1:npoin), vertices_i(1:pnode), point_j(1:ndime), vol_coords_j(1:pnode) )  !< 2016Mar23

      tetracoords_j(ii,1:pnode) = vol_coords_j(1:pnode)
      !-------------------------------------------------------------------||---!
    enddo elementary_loop
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !============================================================================!

  !============================================================================!
  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date
  !> @brief
  !> @details
  !-----------------------------------------------------------------------||---!
  subroutine commdom_plepp_locator_send_nodal_var01(props_i, elemts_i, n_dist_j, tetracoords_j, props_j)
  implicit none
  real(rp),    intent( in) :: props_i(npoin)
  integer(ip), intent( in) :: n_dist_j
  integer(ip), intent( in) ::      elemts_i(n_dist_j)
  real(rp),    intent( in) :: tetracoords_j(n_dist_j,mnode)
  real(rp),    intent(out) ::       props_j(n_dist_j)

  integer(ip) :: ii, ielem, pelty, pnode
  integer(ip) :: vertices_i(mnode)
  real(rp)    :: vol_coords_j(mnode)
  real(rp)    :: prop_j(mnode)

  !-----------------------------------------------------------------------||---!
  if(IMASTER) then
  else
    !---------------------------------------------------------------------||---!
    elementary_loop: &
    do ii = 1,n_dist_j
      ielem = elemts_i(ii)
      pelty = ltype(ielem)
      pnode = nnode(pelty)
      vertices_i(  1:pnode) = lnods(1:pnode,ielem)
      vol_coords_j(1:pnode) = tetracoords_j(ii,1:pnode)
            prop_j(1:pnode) =       props_i( vertices_i(1:pnode) )

      props_j(ii) = dot_product( prop_j(1:pnode), vol_coords_j(1:pnode) )
    enddo elementary_loop
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !============================================================================!

  !============================================================================!
  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date
  !> @brief
  !> @details
  !-----------------------------------------------------------------------||---!
  !subroutine par_commdom_XXX()
  !end subroutine
  !============================================================================!

#endif
end module mod_commdom_plepp
!==============================================================================!
!==============================================================================!

!==============================================================================!
!==============================================================================!
!------------------------------------------------------------------------
!> @addtogroup par_commdom
!> @ingroup    Coupling
!> @{
!>   @name    PLE++ functions
!>   @file    par_commdom_plepp.f90
!>   @author  JM Zavala-Ake
!>   @date
!>   @brief
!>   @details
!> @{
!------------------------------------------------------------------------
!===================================================================================!
!===================================================================================!
