!===================================================================================!
!  I am your father...
! 
!< 2016AUG12 @ Chicago, USA 
!
!===================================================================================!
module mod_mui
  use def_kintyp
  !usa def_parall, only : nproc_par, iproc_par
  use def_master, only : inotmaster, imaster, isequen
  use def_master, only : title, dtime   
  use def_domain, only : coord 
  use def_domain, only : npoin, ndime  
  implicit none 
#ifndef MPI_OFF
  include  'mpif.h'
#endif
  ! 
#ifdef MUI
  external            mui_create_uniface3d_f, mui_push_f, mui_destroy_uniface3d_f
  integer, pointer :: mui_uniface
#endif 
  ! 
  ! Variables
  integer(ip), parameter       :: cp = 16_ip  
  !
  ! Data Structures
  !  
  type     MUI_COUPLING_XX 
    character(cp)              :: name_mesh       = '' 
    integer(ip)                ::   id_mesh       = -1_ip
    !
    integer(ip)                ::    n_props      =  -1_ip 
    integer(ip),      pointer  ::   id_props(:)   => null() 
    integer(ip),      pointer  ::  dim_props(:)   => null()
    integer(ip),      pointer  :: size_props(:)   => null()
    character(cp),    pointer  :: name_props(:)   => null() 
    !
    real(rp),    pointer       ::      props(:,:) => null()
    real(rp),    pointer       ::  positions(:)   => null() 
    integer(ip), pointer       :: vertex_ids(:)   => null() 
    !
    integer(ip)                :: idx_vertices_j  =  -1_ip    
    integer(ip)                ::   n_vertices_j  =  -1_ip      
    integer(ip), pointer       :: idx_coords_j(:) => null()
    !
    real(rp)                   :: tsl             = -huge(1.0_rp) 
    !
    logical(ip)                :: wet_processor   = .false.   
    !
    integer(ip)                :: current_code    =  -1_ip   
    integer(ip)                :: n_modules       =  -1_ip
    integer(ip),      pointer  :: modules(:)      => null()
    !
    logical(ip)                :: debug           =  .true.
    !
  end type MUI_COUPLING_XX  
  ! 
 !type(MUI_COUPLING) PRE_CPLNG 
  !
  type MUI_COUPLING    !< type(COMMDOM_PLEPP_COUPLING) CPLNG
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

    integer(ip), pointer :: idx_coords_j(:)      => null()
    !
    character(cp)  :: namei       = ''
    character(cp)  :: namej       = ''
    character(cp)  :: app_type    = ''
    character(cp)  :: app_name    = ''
    character(cp)  :: module_name = ''
    !
    integer(ip)                ::   current_code  =  -1_ip
    integer(ip)                ::    n_props      =  -1_ip
    integer(ip),      pointer  ::   id_props(:)   => null()
    integer(ip),      pointer  ::  dim_props(:)   => null()
    integer(ip),      pointer  :: size_props(:)   => null()
    integer(ip),      pointer  :: moduls(:)       => null()
    character(cp),    pointer  :: name_props(:)   => null()
    !  
    logical(ip)                :: debug           =  .true.
    ! 
integer(ip) :: code_i, code_j 
integer(ip) ::  module_i, module_j 
integer(ip) ::  what_i, what_j 
integer(ip) ::  fixbo_i, fixbo_j, current_fixbo 
integer(ip) :: n_dof
real(rp)    :: tolerance 
    !
  end type MUI_COUPLING  
  !
  type(MUI_COUPLING), save :: MUI_CPLNG
  !
  private
  !
  public :: MUI_COUPLING, MUI_CPLNG 
  public :: mui_configure_xx 
  public :: mui_initialize
  public :: mui_finalize
  public :: mui_exchange
  public :: mui_init 
  !
contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
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
  !>   call mui_init()      !> JMAKE
  !>   if( PAR_UNIVERSE_SIZE > 1 ) 
  !>   ...
  !>  @endcode
  !> 
  !-----------------------------------------------------------------------||---!
  subroutine mui_init( CPLNG )
  use mod_parall, only: PAR_COMM_CURRENT, PAR_UNIVERSE_SIZE
  use mod_parall, only: PAR_COMM_WORLD, PAR_COMM_UNIVERSE, PAR_COMM_MY_CODE
  implicit none
  type(MUI_COUPLING), intent(inout) :: CPLNG
#ifdef MUI 
  character(cp)  :: send, recv = ''
  character(cp)  :: token
  integer(ip)    :: app_ok, n_send, n_recv
  integer(ip)    :: oki, okj
  integer(4)     :: iarg

  CPLNG%module_name =  ''

  !-----------------------------------------------------------------------||---!
  !app_type = "ALYA_CFD"  !< coupling with syrthes 
  CPLNG%app_type = "MUI" 
  CPLNG%namei    = "DIRIC"     
  CPLNG%namej    = "NEUMA"    
  CPLNG%tol      =  1e-3 
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  call mui_create()

  do iarg = 1, command_argument_count() 
    token = ""
    call get_command_argument(iarg, token)
    call mui_set_argvs(trim(token), len_trim(token))
  enddo

  token = "--name"
  call mui_analyse_argvs(trim(token), len_trim(token))

  CPLNG%app_name = ""
  call mui_get_argvs(CPLNG%app_name)
  CPLNG%module_name = trim(CPLNG%app_name)
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
!!  world_comm = PAR_COMM_UNIVERSE
!!  call commdom_create_commij(world_comm, CPLNG%local_comm)
  CPLNG%local_comm = MPI_COMM_NULL  
  call mui_create_uniface3d_f(trim(CPLNG%app_name), CPLNG%local_comm)

!! 
!! SEE: uniface.h
!!        std::unique_ptr<communicator> comm;
!!      comm_mpi.h
!!        MPI_Comm domain_local_;
!!
!!
  PAR_COMM_WORLD    = CPLNG%local_comm
  PAR_COMM_CURRENT  = CPLNG%local_comm
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  CPLNG%commij  = -1
  CPLNG%commji  = -1
  PAR_UNIVERSE_SIZE = -1 !=> if( PAR_UNIVERSE_SIZE > 1 ) ... 
  !-----------------------------------------------------------------------||---!
#endif 
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine mui_finalize( CPLNG )
  implicit none
  type(MUI_COUPLING),       intent(inout)   :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef MUI
  call mui_destroy_uniface3d_f()
  call mui_delete()
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine mui_configure_xx( CPLNG ) 
  use def_master, only : current_code
  use def_master, only : ID_TEMPER
  implicit none
  type(MUI_COUPLING),       intent(inout)   :: CPLNG
  !
  character(50)   ::  corsc
  integer(ip)     ::  kfl_required ! 1: yes, 0: no
  !
  corsc(:) = ' '
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef MUI
    !
    CPLNG%n_props        = 1  
    !
    if( .not.associated(CPLNG%id_props)    ) allocate( CPLNG%id_props(CPLNG%n_props)   ) 
    if( .not.associated(CPLNG%dim_props)   ) allocate( CPLNG%dim_props(CPLNG%n_props)  )
    if( .not.associated(CPLNG%size_props)  ) allocate( CPLNG%size_props(CPLNG%n_props) )
    if( .not.associated(CPLNG%name_props)  ) allocate( CPLNG%name_props(CPLNG%n_props) ) 
    !
    if(trim(title) == "DIRICH") then 
      CPLNG%current_code   =  1
      !
      CPLNG%id_props(1)    = -1_ip  
      CPLNG%dim_props(1)   =  1_ip  
      CPLNG%name_props(1)  = 'TEMPE'
      !
      CPLNG%moduls(1)      =  ID_TEMPER 
    endif 
    !
    if(trim(title) == "NEUMAN") then  
      CPLNG%current_code   =  2
      !
      CPLNG%id_props(1)    = -1_ip
      CPLNG%dim_props(1)   =  1_ip 
      CPLNG%name_props(1)  = 'TFLUX'
      !
      CPLNG%moduls(1)      =  ID_TEMPER 
    endif 
    !
    CPLNG%current_code     = current_code 
    !
    if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then
      print *, "[mui_init] ", "'", trim(title), "."
    endif
    !
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine mui_initialize( CPLNG )
  use def_domain, only : ltype, lnods, coord, nelem 
  implicit none
  type(MUI_COUPLING),       intent(inout)   :: CPLNG
  !
  integer(ip)          :: id 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef MUI
    !
    !
    if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then
    endif
    !
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine mui_exchange( CPLNG, nameij, nameji, tij, tji, varij, varji ) 
  use def_domain,           only: coord, mnode, ndime, npoin
  implicit none
  type(MUI_COUPLING),              intent(inout)   :: CPLNG
  character(4), optional,          intent(in   )   :: nameij  
  character(4), optional,          intent(in   )   :: nameji 
  real(rp),     optional,          intent(in   )   :: tij  
  real(rp),     optional,          intent(in   )   :: tji 
  !
  real(rp),     optional,          intent(in   )   :: varij(npoin)
  real(rp),     optional,          intent(out  )   :: varji(npoin)
  ! 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip)   :: nij=0.0, nji=0.0, k
  real(rp)      :: rij(3)
  real(rp)      :: rji(3)

if( INOTMASTER ) nij=npoin
if( INOTMASTER ) nji=npoin

#ifdef MUI
  !-----------------------------------------------------------------------||---!
  ! push data to the other solver
  if( present(varij) ) then 
    print *, "[mui_exchange] ",  "'"//trim(title)//"."//trim(nameij)//"'", " -->"
    do k = 1,nij 
      rij(1:3)     = 0.0 
      rij(1:ndime) = coord(1:ndime,k)
      call mui_push_f( nameij, rij(1), rij(2), rij(3), varij(k) )
    enddo 
    call mui_commit_f( tij )
  endif 
  !-----------------------------------------------------------------------||---!
  ! fetch data from the other solver
  if( present(varji) ) then
    print *, "[mui_exchange] ",  "'"//trim(title)//"."//trim(nameji)//"'", "<-- "
    do k = 1,nji 
       call mui_fetch_moving_average_mean_f(                 & 
                                                   nameji,       & 
                                                   rji(1), rji(2), rji(3), tji, & 
                                                   varji(k) &  
                                                 )   
    enddo 
  endif
  !-----------------------------------------------------------------------||---!
  if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then  
!    print *, "[mui_exchange] ", "'"//trim(title)//"'", "'", trim( CPLNG%name_props(id) ), "'"   
  endif
  !-----------------------------------------------------------------------||---!
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine mui_set_mesh( id_fixbo_j, stride )
  use def_domain, only: coord, mnode, nelem, ndime, npoin
  use def_domain, only: ltype, nnode, ngaus, lnods, coord
  implicit none
  integer(ip), optional, intent(in) :: id_fixbo_j
  integer(ip), optional, intent(in) :: stride
  ! 
  integer(ip)   :: n_dime_i, n_node_i
  integer(ip)   :: n_vertices_i, n_elements_i, n_vertices_j
  !
  real(rp),    pointer :: vertex_coords_i(:,:) => null()
  integer(ip), pointer ::    vertex_num_i(:,:) => null()
  integer(ip), pointer ::   vertex_type_i(:  ) => null()
  real(rp),    pointer :: vertex_coords_j(:,:) => null()
  integer(ip), pointer :: interior_list_j(:  ) => null()
  !-----------------------------------------------------------------------||---!
  if(.True.) then
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
!      if(id_fixbo_j>0) call mui_on_field(PLEPP_CPLNG, id_fixbo_j, n_vertices_j, vertex_coords_j)
    endif
    !---------------------------------------------------------------------||---!

    !---------------------------------------------------------------------||---!
    !
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| PRIVATE |---!
!-------------------------------------------------------------------------||---!


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
  subroutine mui_on_field(CPLNG, where_number)
  use def_domain, only:  nfiel, kfl_field, xfiel
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  type(MUI_COUPLING),       intent(inout)   :: CPLNG
  integer(ip),                  intent(in   )   :: where_number
  !
  integer(ip) :: ipoin, ifiel, idime, i_coords_j, dummy 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(nfiel==0) then
    print *, "[commdom_plepp_on_field]", " Into '*.dom.dat' include -> 'FIELD=1, DIMENSION=1, NODE'"
    call runend('EXIT!!')
  endif
  !
  if(where_number>0) then 
    ifiel = 1_ip
    idime = 1_ip
    CPLNG%n_vertices_j = 0_ip
    !
    if(INOTMASTER) then
      !
      do ipoin =1_ip,npoin
        if( int(xfiel(ifiel) % a(idime,ipoin,1)) == where_number) then
        !!print*, coord(:,ipoin), xfiel(ifiel) % a(idime,ipoin)  
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
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
!    if( .not.associated(CPLNG%idx_coords_j) ) allocate( CPLNG%idx_coords_j(CPLNG%n_vertices_j) )
!    if( .not.associated(CPLNG%positions)    ) allocate( CPLNG%positions(ndime*CPLNG%n_vertices_j) )
!    CPLNG%idx_coords_j = 0.0_rp
!    CPLNG%positions    = 0.0_rp
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER) then
      dummy      = 0_ip 
      i_coords_j = 0_ip
      do ipoin = 1,npoin
        if( int(xfiel(ifiel) % a(idime,ipoin,1)) == where_number) then
         i_coords_j = i_coords_j + 1_ip
         !
!         CPLNG%idx_coords_j(i_coords_j) = ipoin
         do idime = 1,ndime 
           dummy  = dummy + 1
!           CPLNG%positions(dummy) = coord(idime,ipoin)  
         enddo 
         !
       endif
      enddo
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---! 
  if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) print *, "[mui_on_field]", " n_coords_j:", CPLNG%n_vertices_j, "where:", where_number
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine mui_xxx( CPLNG )
  implicit none
  type(MUI_COUPLING),       intent(inout)   :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef MUI
  if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then 
    print *, "[mui_xxx]"
  endif 
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
end module mod_mui
