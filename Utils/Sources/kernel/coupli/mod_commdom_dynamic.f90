!==============================================================================!
  !
  !< 2015Mar20 -> created (from 'mod_commdom_plepp')
  !< 2015Jul03 -> add 'commdom_dynamic_reduce_sum'
  !< 2015AGO13 -> add 'commdom_dynamic_kdtree_01'
  !< 2017JAN07
  !< 2017JAN09
  !
  !-----------------------------------------------------------------------||---!
  !
  !   + current_code                                      ___________current_task
  !   |_Alya                                       ______|_____
  !     |_call Turnon()                            ITASK_TURNON 02
  !     |_call Iniunk()                            ITASK_INIUNK 03
  !     |_time: do while
  !       |_call Timste()                          ITASK_TIMSTE 04
  !       |_reset: do
  !         |_call Begste()                        ITASK_BEGSTE 05
  !           |_block: do while
  !             |_coupling: do while
  !               |_call Begzon()                  ITASK_BEGZON 19
  !               |_modules: do while                               / TASK_BEGITE  14
  !                 |_call Doiter()                ITASK_DOITER 06-|
  !                 |_call Concou()                ITASK_CONCOU 07  \_ITASK_ENDITE 15
  !               |_call Endzon()                  ITASK_ENDZON 20
  !             |_call Conblk()                    ITASK_CONBLK 08
  !       |_call Endste()                          ITASK_ENDSTE 10
  !
  !-----------------------------------------------------------------------||---!
  !          __
  ! BLOCK 3_   |
  !   1 X   |  |--current_block
  !   2 Y Z |  |
  !   3 W  _|-----current_module
  ! END_BLOCK__|
  !
  !-----------------------------------------------------------------------||---!
  !
  ! <code, block, modul, task, when, send|recv>
  !
  ! modules: ID_KERNEL=0, ID_NASTIN=1, ID_TEMPER=2, ID_NASTAL=6, ID_ALEFOR=7, ID_SOLIDZ=10
  !    when: ITASK_BEFORE=1, ITASK_AFTER=2
  !
  !-----------------------------------------------------------------------||---!
!==============================================================================!
module mod_commdom_dynamic
  use mod_commdom_alya,  only: INONE
#ifdef COMMDOM
  use def_parame,        only: ip, rp
  use def_master,        only: inotmaster, isequen
  use def_domain,        only: coord, mnode, nelem, ndime, npoin
  use def_domain,        only: LESET
  use def_domain,        only: ltype, nnode, ngaus, lnods, coord
  use mod_commdom_plepp, only: COMMDOM_PLEPP_COUPLING
  use mod_std
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  logical(ip),   parameter :: STATISTICS = .False.                                       !<- ON|OFF. 2017JAN07
  !
  integer(ip),   parameter :: cp   = 64
  character(cp), parameter :: frmt = '(E11.4)'
  !
  logical(ip) :: initiated = .false.
  logical(ip) :: debug     = .true.
  !
  type COMMDOM_DYNAMIC_FIXNO
    integer(ip)          :: n_fixno        =  0
    integer(ip)          :: n_dof          =  0
    logical(ip)          :: initiated      = .false.
    logical(ip), pointer ::  launched(:  ) => null()
    integer(ip), pointer ::     fixno(:,:) => null()
  end type COMMDOM_DYNAMIC_FIXNO
  !
  type(COMMDOM_DYNAMIC_FIXNO) FIXNO
  !
  !-----------------------------------------------------------------------||---!
  type COMMDOM_DYNAMIC_RELAXATION
    real(rp),    pointer :: var_ij(:,:) => null()
  end type COMMDOM_DYNAMIC_RELAXATION
  !
  type(COMMDOM_DYNAMIC_RELAXATION) RELAXATION
  !
  !-----------------------------------------------------------------------||---!
  type COMMDOM_VERTEX_PROPERTIES
    real(rp),    pointer :: vertex_props_j(:,:) => null()
    integer(ip)          :: ndime_props         = 0_ip
    integer(ip)          :: n_pts               = 0_ip
  end type COMMDOM_VERTEX_PROPERTIES
  !
  type(COMMDOM_VERTEX_PROPERTIES) CPLNG_PROPS
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  private
    public :: CPLNG_PROPS
    public :: commdom_dynamic_deallocate
    public :: commdom_dynamic_set_mesh
    public :: commdom_dynamic_exchange02
    public :: commdom_dynamic_outvar
    public :: commdom_dynamic_check_fixno
    public :: commdom_dynamic_set_vals
    public :: commdom_dynamic_reduce_sum
    public :: commdom_dynamic_kdtree
    public :: STATISTICS
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_set_mesh( id_fixbo_j, stride )
  use mod_commdom_plepp,  only: PLEPP_CPLNG
  use mod_commdom_aitken, only: RELAXATION
  use mod_commdom_aitken, only: commdom_aitken_init, commdom_aitken_allocate
  implicit none
  integer(ip), optional, intent(in) :: id_fixbo_j
  integer(ip), optional, intent(in) :: stride
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( .not.initiated ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    call commdom_dynamic_create(                  PLEPP_CPLNG, id_fixbo_j )
    call commdom_dynamic_allocate(                PLEPP_CPLNG, stride     )
    call commdom_dynamic_get_coupling_properties( PLEPP_CPLNG, id_fixbo_j )

    call commdom_dynamic_statistics( PLEPP_CPLNG%n_ij, PLEPP_CPLNG%n_ji   ) !< 2016DIC21
!
!if(n_send>0) call commdom_dynamic_kdtree_01( n_send, PLEPP%dist_coords_j ) !<  2015AGO13
!
    call commdom_dynamic_fixno_allocate(          PLEPP_CPLNG, stride     )
    !
    initiated = .true.
    !
    !print*, "[commdom_dynamic_set_mesh]"
    !
   !call commdom_aitken_init(     RELAXATION, PLEPP_CPLNG%n_ji, PLEPP_CPLNG%stride, 3_ip, 0.8_rp)
   !call commdom_aitken_allocate( RELAXATION )
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_create( PLEPP, id_fixbo_j )
  use mod_commdom_plepp, only: commdom_plepp_on_field
  use def_master,        only: displ
  use mod_communications, only :  PAR_MAX
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  integer(ip), optional, intent(in) :: id_fixbo_j
  !
  real(rp),    pointer :: vertex_coords_i(:,:) => null()
  integer(ip), pointer ::    vertex_num_i(:,:) => null()
  integer(ip), pointer ::   vertex_type_i(:  ) => null()
  real(rp),    pointer :: vertex_coords_j(:,:) => null()
  integer(ip), pointer :: interior_list_j(:  ) => null()
  !
  integer(ip) :: n_vertices_i, n_elements_i, n_vertices_j
  integer(ip) :: n_dime_i, n_node_i
  !
  integer(ip) :: element_type
  !
  n_vertices_i = 0
  n_elements_i = 0
  n_vertices_j = 0
  n_dime_i     = 0
  n_node_i     = 0
  !
  element_type = -1
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !----------------------------------------------------------| init_mesh |---!
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
      CPLNG_PROPS%vertex_props_j => displ(:,:,1)     !< 2017Nov19
      CPLNG_PROPS%ndime_props    = ndime             !< 2017Nov19

    else
      allocate( vertex_coords_i(0,0) )
      allocate(    vertex_num_i(0,0) )
      allocate(   vertex_type_i(0  ) )
      allocate( vertex_coords_j(0,0) )
    endif
    !
    if( present(id_fixbo_j) ) then
      if(id_fixbo_j>0) call commdom_plepp_on_field(PLEPP,    id_fixbo_j, &
                                                           n_vertices_j, &
                                                           vertex_coords_j)
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !-------------------------------------------------------| ELEMENT TYPE |---!
    if(all(vertex_type_i==30).or.(all(vertex_type_i==10))) element_type = 4 ! TETRAS=30 | TRIA=10
    if(all(vertex_type_i==37)                            ) element_type = 6 ! HEXAS=37
    call PAR_MAX(element_type,'IN MY CODE')
    !print *, trim(PLEPP%module_name), " element_type:", element_type

    !-----------------------------------------------------| create_locator |---!
    if( .not.(PLEPP%commij == -1) ) then
      call commdom_locator_create2(PLEPP%local_comm, PLEPP%commij, PLEPP%tol, element_type)
    endif
    if( .not.(PLEPP%commji == -1) ) then
      call commdom_locator_create2(PLEPP%local_comm, PLEPP%commji, PLEPP%tol, element_type)
    endif
    !
    call commdom_locator_set_cs_mesh(   n_vertices_i, &
                                        n_elements_i, &
                                     vertex_coords_i, &
                                        vertex_num_i, &
                                       vertex_type_i, &
                                        n_vertices_j, &
                                     vertex_coords_j, &
                                               ndime, &
                          CPLNG_PROPS%vertex_props_j, &
                             CPLNG_PROPS%ndime_props )  !< 2016Ago25
    !
    !call commdom_locator_save_dist_coords(0, PLEPP%local_comm)
    !
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_deallocate( PLEPP )
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( ( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ).and.initiated ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    !
    call commdom_locator_destroy()
    !
    if( associated( PLEPP%tetra_coords_j )  ) deallocate( PLEPP%tetra_coords_j   )
    if( associated( PLEPP%dist_locations_i) ) deallocate( PLEPP%dist_locations_i )
    if( associated( PLEPP%dist_coords_j)    ) deallocate( PLEPP%dist_coords_j    )
    if( associated( PLEPP%var_ij)           ) deallocate( PLEPP%var_ij           )
    if( associated( PLEPP%var_ji)           ) deallocate( PLEPP%var_ji           )
    if( associated( PLEPP%interior_list_j)  ) deallocate( PLEPP%interior_list_j  )
    if(CPLNG_PROPS%ndime_props > 0_ip ) then
      if( associated( PLEPP%dist_props_j)   ) deallocate( PLEPP%dist_props_j     )
    end if
    !
    PLEPP%n_ij   = 0
    PLEPP%n_ji   = 0
    !
    initiated = .false.
    !
    !print*, "[commdom_dynamic_deallocate]
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_allocate( PLEPP, stride )
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  integer(ip),        optional, intent(in   ) :: stride
  !
  integer(ip)   :: n_recv, n_send
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    n_recv = 0
    n_send = 0
    call commdom_locator_get_n_dist_points( n_send )
    call commdom_locator_get_n_interior(    n_recv )
    !
    if( .not.associated( PLEPP%tetra_coords_j )  ) allocate( PLEPP%tetra_coords_j(            n_send, mnode) )
    if( .not.associated( PLEPP%dist_locations_i) ) allocate( PLEPP%dist_locations_i(          n_send       ) )
    if( .not.associated( PLEPP%dist_coords_j)    ) allocate( PLEPP%dist_coords_j(             n_send*ndime ) )
    if( .not.associated( PLEPP%var_ij)           ) allocate( PLEPP%var_ij(            stride, n_send       ) )
    if( .not.associated( PLEPP%var_ji)           ) allocate( PLEPP%var_ji(            stride, n_recv       ) )
    if( .not.associated( PLEPP%interior_list_j)  ) allocate( PLEPP%interior_list_j(           n_recv       ) )
    if(CPLNG_PROPS%ndime_props > 0_ip ) then
      if( .not.associated( PLEPP%dist_props_j)   ) allocate( PLEPP%dist_props_j( n_send * CPLNG_PROPS%ndime_props  ) )
    end if
    !
    PLEPP%n_ij   = n_send
    PLEPP%n_ji   = n_recv
    PLEPP%stride = stride
    !
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_get_coupling_properties( PLEPP, id_fixbo_j )
  use mod_commdom_plepp, only: commdom_plepp_locator_send_nodal_var00
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  integer(ip), pointer :: interior_list_j(:  ) => null()
  integer(ip), optional, intent(in) :: id_fixbo_j
  integer(ip)   :: n_recv, n_send
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    n_recv = 0
    n_send = 0
    call commdom_locator_get_n_dist_points( n_send )
    call commdom_locator_get_n_interior(    n_recv )
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER) then
      PLEPP%dist_locations_i = -1
      PLEPP%dist_coords_j    = -1
      call commdom_locator_get_dist_locations( PLEPP%dist_locations_i )
      call commdom_locator_get_dist_coords(    PLEPP%dist_coords_j    )
      call commdom_plepp_locator_send_nodal_var00( PLEPP%dist_locations_i, &
                                                   PLEPP%dist_coords_j,    &
                                                   n_send, &
                                                   PLEPP%tetra_coords_j )
      if(CPLNG_PROPS%ndime_props > 0_ip ) then
        call commdom_locator_get_dist_props(    PLEPP%dist_props_j    )
      end if
      !
      if( present(id_fixbo_j) ) then
        if( id_fixbo_j>0)  then
          allocate( interior_list_j( n_recv  ) )
          interior_list_j(1:n_recv) = -1
          call commdom_locator_get_interior_list( interior_list_j(1:n_recv) )
          PLEPP%interior_list_j(1:n_recv) = PLEPP%idx_coords_j( interior_list_j(1:n_recv)  )
          deallocate( interior_list_j )
        else
          PLEPP%interior_list_j(1:n_recv) = -1
          call commdom_locator_get_interior_list( PLEPP%interior_list_j(1:n_recv) )
        endif
      endif
      !
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!---------------------------------------------------------------| EXCHANGE |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_exchange02(prop_i, prop_j, stride, debug)
  use mod_commdom_plepp, only: PLEPP_CPLNG
  implicit none
  real(rp),              intent( in) :: prop_i(:,:)
  real(rp),              intent(out) :: prop_j(:,:)
  integer(ip),           intent( in) :: stride
  logical(ip),           intent( in) :: debug
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_dynamic_exchange02_00(PLEPP_CPLNG, prop_i, prop_j, stride, debug)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_exchange02_00(PLEPP, prop_i, prop_j, stride, debug)
  use mod_commdom_plepp, only: commdom_plepp_locator_send_nodal_var01
  implicit none
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  real(rp),                     intent(in   ) :: prop_i(:,:)
  real(rp),                     intent(out  ) :: prop_j(:,:)
  integer(ip),                  intent(in   ) :: stride
  logical(ip),                  intent(in   ) :: debug
  !
  real(rp)      :: daux(4)
  character(cp) :: saux(4)
  integer(ip)   :: ii
!
!real(rp)      :: prop_in, prop_out
!prop_out = huge(1.0_rp)
!prop_in  = huge(1.0_rp)
!
  if( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if((stride<1_ip).or.(stride>PLEPP%stride)) then
      print *, "[commdom_plepp_exchange02] ERROR: stride<PLEPP%stride:", stride, "<",PLEPP%stride
      call runend('EXIT!!')
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    PLEPP%var_ij(1_ip:stride,1_ip:PLEPP%n_ij) = -1
    !
    if(INOTMASTER) then
      do ii = 1,stride
        call commdom_plepp_locator_send_nodal_var01(prop_i(ii,1:npoin),                                                 & !<---
                                                    PLEPP%dist_locations_i(  1_ip:PLEPP%n_ij           ), PLEPP%n_ij, &
                                                    PLEPP%tetra_coords_j(    1_ip:PLEPP%n_ij,1_ip:mnode),               & !--->
                                                    PLEPP%var_ij(         ii,1_ip:PLEPP%n_ij           ) )
        !-----------------------------------------------------------------||---!
        !                                                                      !
        !-----------------------------------------------------------------||---!
        daux(1) = minval( PLEPP%var_ij(ii,1_ip:PLEPP%n_ij) )
        if( PLEPP%n_ij /= 0_ip ) then
           daux(2) = sum( PLEPP%var_ij(ii,1_ip:PLEPP%n_ij) )/real(PLEPP%n_ij,rp)
        else
           daux(2) = 0.0_rp
        end if
        daux(3) = maxval( PLEPP%var_ij(ii,1_ip:PLEPP%n_ij) )
        !
        write(saux(1), frmt) daux(1)
        write(saux(2), frmt) daux(2)
        write(saux(3), frmt) daux(3)
        write(saux(4), '(I2)') ii
        !
        if(PLEPP%n_ij/=0 .and. debug) &
          print*,"["//trim(PLEPP%module_name)//trim(saux(4))//"]->["//trim(saux(1))//","//trim(saux(2))//","//trim(saux(3))//"]"
      enddo
!      prop_in = daux(1)
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
! 2015Jul03
!
!    if(.not.(PLEPP%commij==-1)) then
!      call commdom_reduce_min_real(prop_in, prop_out, PLEPP%local_comm, PLEPP%commij)
!    endif
!    if(.not.(PLEPP%commji==-1)) then
!      call commdom_reduce_min_real(prop_in, prop_out, PLEPP%local_comm, PLEPP%commji)
!    endif
!
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    PLEPP%var_ji(1_ip:stride,1_ip:PLEPP%n_ji) = -1
    call commdom_locator_exchange_double_stride(PLEPP%var_ij, PLEPP%var_ji, stride) !> send, recv
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER) then
    do ii = 1,stride
      prop_j(ii,PLEPP%interior_list_j) = PLEPP%var_ji(ii,1:PLEPP%n_ji)
      !
      daux(1) = minval( prop_j(ii,PLEPP%interior_list_j) )
      if( PLEPP%n_ji /= 0_ip ) then
         daux(2) = sum( prop_j(ii,PLEPP%interior_list_j) )/real(PLEPP%n_ji,rp)
      else
         daux(2) = 0.0_rp
      end if
      daux(3) = maxval( prop_j(ii,PLEPP%interior_list_j) )
      !
      write(saux(1), frmt) daux(1)
      write(saux(2), frmt) daux(2)
      write(saux(3), frmt) daux(3)
      write(saux(4), '(I2)') ii
      !
      if(PLEPP%n_ji/=0 .and. debug ) &
        print *, "["//trim(PLEPP%module_name)//trim(saux(4))//"]<-["//trim(saux(1))//","//trim(saux(2))//","//trim(saux(3))//"]"
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
  subroutine commdom_dynamic_reduce_sum(prop_in, prop_out, debug)              !< 2015Jul03
  use mod_commdom_plepp, only: PLEPP_CPLNG
  implicit none
  real(rp),              intent(in ) :: prop_in
  real(rp),              intent(out) :: prop_out
  logical(ip), optional, intent(in ) :: debug
  !
  prop_out = huge(1.0_rp)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( .not.(PLEPP_CPLNG%commij==-1) ) call commdom_reduce_sum_real(prop_in, prop_out, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commij)
  if( .not.(PLEPP_CPLNG%commji==-1) ) call commdom_reduce_sum_real(prop_in, prop_out, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commji)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_set_vals(PHI_1, PHI_2, relax_op, res2, debug)
  use mod_commdom_plepp, only: PLEPP_CPLNG
  implicit none
  real(rp),              intent(in ) :: PHI_1(npoin) ! PHI_{m  }^{n}
  real(rp),              intent(out) :: PHI_2(npoin) ! PHI_{m+1}^{n}
  real(rp),    optional, intent(in ) :: relax_op
  real(rp),    optional, intent(out) :: res2(2)
  logical(ip), optional, intent(in ) :: debug
  !
  real(rp) :: relax = 1.0_rp
  if(present(relax_op)) relax = relax_op
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_dynamic_set_vals_00(PLEPP_CPLNG, PHI_1, PHI_2, relax, res2, debug)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_set_vals_00(PLEPP, PHI_1, PHI_2, relax, modR2, debug)
  use def_master,          only: ittim
  implicit none
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  ! PHI_{m}^{n}; n:time loop, m:iteration loop
  !
  real(rp),              intent(in ) :: PHI_1(npoin) ! PHI_{m  }^{n}
  real(rp),              intent(out) :: PHI_2(npoin) ! PHI_{m+1}^{n}
  real(rp),              intent(in ) :: relax
  real(rp),    optional, intent(out) :: modR2(2)
  logical(ip), optional, intent(in ) :: debug
  !
  !real(rp) :: modR2 = 0.0_rp
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER) then
      !
      !! PHI_{m+1}^{n} = omega PHI_{m+1}^{n} + (1-omega) PHI_{m}^{n}
      !PHI_2(PLEPP%interior_list_j) = PHI_2(PLEPP%interior_list_j) * (    relax) + &
      !                               PHI_1(PLEPP%interior_list_j) * (1.0-relax)
      PHI_2(PLEPP%interior_list_j) = &
                                     PHI_2(PLEPP%interior_list_j) + &
                           relax * ( PHI_1(PLEPP%interior_list_j) - PHI_2(PLEPP%interior_list_j) )
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    !  R_{m+1}^{n} = PHI_{m+1}^{n} - PHI_{m}^{n}
    if( present(modR2) ) then
      modR2 = 0.0
      if(INOTMASTER) then
        modR2(2) = dot_product( PHI_2(PLEPP%interior_list_j) - PHI_1(PLEPP%interior_list_j) ,&
                                PHI_2(PLEPP%interior_list_j) - PHI_1(PLEPP%interior_list_j) )
        modR2(1) = dot_product( PHI_2(PLEPP%interior_list_j),  PHI_2(PLEPP%interior_list_j) )
      endif
      if(ittim<1) modR2 = huge(1.0_rp)
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER) then
      if( present(debug).and.debug ) then
        print *, "[commdom_dynamic_set_vals]", minval( PHI_1( PLEPP%interior_list_j) ), &
                                                  sum( PHI_1( PLEPP%interior_list_j) )/PLEPP%n_ji, &
                                               maxval( PHI_1( PLEPP%interior_list_j) ), "<--"
        if(relax/=1.0) &
        print *, "[commdom_dynamic_set_vals]", minval( PHI_2(PLEPP%interior_list_j) ), &
                                                  sum( PHI_2(PLEPP%interior_list_j) )/PLEPP%n_ji, &
                                               maxval( PHI_2(PLEPP%interior_list_j) ), "-->"
      endif
      !-------------------------------------------------------------------||---!
      !                                                                        !
      !-------------------------------------------------------------------||---!
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| ++++ |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_check_fixno(fixnode, idofn, fixval, ToDo)
  use mod_commdom_plepp, only: PLEPP_CPLNG
  implicit none
  integer(ip),  intent(inout), pointer  :: fixnode(:,:)
  integer(ip),  intent(in   )           :: idofn
  integer(ip),  intent(in   )           :: fixval
  logical(ip),  intent(in   )           :: ToDo
  !
  !integer(ip) :: idofn
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  !do idofn = 1,FIXNO%n_dof
    !
    if( FIXNO%n_dof<idofn.or.idofn<=0 ) then
      print *, "[commdom_dynamic_check_fixno] ERROR: FIXNO%n_dof/=idofn:",  FIXNO%n_dof, "/=", idofn
      call runend('EXIT!!')
    endif
    !
    if( FIXNO%launched(idofn) ) then
      ! (1) reset fixno
      fixnode(    idofn,1:npoin) = FIXNO%fixno(idofn,1:npoin)
    else
      ! (0) save original fixno (wich dont evolve!!)
      FIXNO%fixno(idofn,1:npoin) = fixnode(    idofn,1:npoin)
      FIXNO%launched(idofn) = .true.
    endif
    !
    ! (2) original fixno + new nodes located
    call commdom_dynamic_check_fixno_00(   PLEPP_CPLNG, fixnode, idofn, fixval, ToDo ) ! fixno, idofn, fixval, ToDo
  !enddo
  !
  ! (3) be happy!!
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_check_fixno_00(PLEPP, fixno, &
                                            idofn, fixval, ToDo)
  use def_domain, only: kfl_codno
!
!use def_kintyp,           only: soltyp
!use def_master,           only: momod, modul
!
  implicit none
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  integer(ip),  intent(inout), pointer  :: fixno(:,:)
  integer(ip),  intent(in   )           :: idofn
  integer(ip),  intent(in   )           :: fixval
  logical(ip),  intent(in   )           :: ToDo
  !
  integer(ip) :: n_fixval, n_coupling
  integer(ip) :: icplng, ipoin, ii
  !
!integer(ip),  pointer :: fixno(:,:)
!type(soltyp), pointer :: solve(:)
!solve                => momod(modul) % solve(1:)
!solve(1) % kfl_fixno => fixno
  !
  if( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ) then
    !
    n_fixval   = 0
    n_coupling = PLEPP%n_ji
    !
    if(INOTMASTER) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    !
    ! n_fixval: cuantos 'fixno' iguales a 'fixval' existen?
    !
    n_fixval = count( fixno(idofn,PLEPP%interior_list_j) == fixval , KIND=ip)
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
        !
        do icplng = 1,n_coupling
          ipoin = PLEPP%interior_list_j(icplng)
          if( fixno(idofn,ipoin) /= fixval ) then
            ii = ii+ 1
           !print *, ii, ",", ipoin, ",<", fixno(:,ipoin), ">, <", kfl_codno(:,ipoin), ">",  fixno(idofn,ipoin),"->", fixval
            fixno(idofn,ipoin) = fixval
          endif
        enddo
        !
        !if(PLEPP%commij /= -1) then
        !  print *, "[commdom_dynamic_check_fixno_i] ", "[i, ipoin, FIXNO, <CODNO>]", idofn
        !  print *, "[commdom_dynamic_check_fixno_i] ","n_coupling-n_fixno:", n_coupling-n_fixval, ", USE FIXNO->", fixval
        !else if(PLEPP%commji /= -1) then
        !  print *, "[commdom_dynamic_check_fixno_j] ", "[j, ipoin, FIXNO, <CODNO>]", idofn
        !  print *, "[commdom_dynamic_check_fixno_j] ","n_coupling-n_fixno:", n_coupling-n_fixval, ", USE FIXNO->", fixval
        !endif
       !!call runend('EXIT!!')
        !
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
  subroutine commdom_dynamic_fixno_allocate( PLEPP, n_dof)
  implicit none
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  integer(ip)                 , intent(in   ) :: n_dof
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ) then
  !
  if(.not.FIXNO%initiated) then
    FIXNO%n_fixno = 0
    FIXNO%n_dof   = n_dof
    if(INOTMASTER) FIXNO%n_fixno = npoin
    if( .not.associated(FIXNO%fixno   )  ) allocate( FIXNO%fixno(   n_dof,npoin) )
    FIXNO%fixno = 0
    !
    if( .not.associated(FIXNO%launched)  ) allocate( FIXNO%launched(n_dof      ) )
    FIXNO%launched  = .false.
    !
    FIXNO%initiated = .true.
    !
  endif
  !
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_set_source_nodes(PLEPP, prop_out)
  implicit none
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  logical(ip), intent(out)  :: prop_out(npoin)
  !
  if( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(INOTMASTER) then
      prop_out(1:npoin              ) = .false.
      prop_out(PLEPP%interior_list_j) = .true.
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  !
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| AUXs |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  !      postp(1) % wopos (1,--) = 'PLEPP'
  !      postp(1) % wopos (2,--) = 'SCALA'
  !
  !----------------------------------------------------------------------------!
  subroutine commdom_dynamic_outvar( )
  use def_master, only: gesca
  use mod_commdom_plepp, only: commdom_plepp_set_source_nodes
  implicit none
  logical(ip), pointer :: touched(:)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( INOTMASTER ) then
      allocate( touched(npoin) )
      !
      call memgen(0_ip, npoin, 0_ip)
      call commdom_plepp_set_source_nodes( touched )
      !
      gesca(1:npoin) = -1.0_rp
      where( touched )  gesca(1:npoin) = 1.0_rp
      !
      deallocate( touched )
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| KDTREE |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_kdtree( prop )
  use mod_commdom_plepp,  only: PLEPP_CPLNG
  implicit none
  real(rp),    intent(inout)        :: prop(:,:)
  !
  real(rp),    pointer              :: dist_coords_j(:) => null()
  real(rp),    pointer              ::   coord_j(:,:,:) => null()
  real(rp),    pointer              ::      R_perp(:,:) => null()
  !
  integer(ip) :: ipoin, idime, stride, n_send, n_recv
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  n_send        =  PLEPP_CPLNG%n_ij
  n_recv        =  PLEPP_CPLNG%n_ji
  stride        =  PLEPP_CPLNG%stride
  dist_coords_j => PLEPP_CPLNG%dist_coords_j
  !
  if(stride /= ndime) call runend("jaja..")
  !
  allocate( coord_j(ndime,n_send,2) )
  allocate(  R_perp(ndime,n_recv  ) )
  !
  if(inotmaster.and.(n_send>0)) then
    !
    ! <- (1) D=> r_gamma => N
    do ipoin = 1,n_send
      do idime = 1,ndime
        coord_j(idime,ipoin,1) = dist_coords_j( ndime*(ipoin-1)+idime )
      enddo
    enddo
    !
    ! -> (2) N: r_gamma -> r_perp
    call commdom_dynamic_kdtree_00(                         n_send, &
                                       coord_j(1:ndime,1:n_send,1), &
                                    Pn=coord_j(1:ndime,1:n_send,2)  )
    !
  endif
  !
  ! (3) N => r_perp => D
  call commdom_locator_exchange_double_stride( coord_j(1:ndime,1:n_send,2), &
                                                R_perp(1:ndime,1:n_recv  ),  ndime )
  !
  if(inotmaster.and.(n_recv>0)) then
    prop(1:ndime,PLEPP_CPLNG%interior_list_j) = R_perp(1:ndime,1:n_recv) - coord(1:ndime,PLEPP_CPLNG%interior_list_j)
    !
    ! (4) D => K u = f
    ! (5) D => sigma/res => N
    ! (6) N  => K u = f
    ! (1) o
    !
  endif
  !
  deallocate( coord_j )
  deallocate( R_perp  )
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_kdtree_00(n_send, coords_j, Dn, Pn )
  use mod_kdtree,        only: kdtree, dpopar
  use def_domain,        only: mnodb
  use def_master,        only: netyp, title
  use def_domain,        only: npoib,nboun,lnodb,ltypb
  use def_master,        only: displ
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip), intent(in)            :: n_send
  real(rp),    intent(in)            :: coords_j(:,:) !coords_j(ndime, npoin)
  real(rp),    intent(out), optional :: Dn(  :)
  real(rp),    intent(out), optional :: Pn(:,:)

  integer(ip)                       :: iboun, inodb, ipoin, idime
  integer(ip)                       :: pblty, pnodb

  !definitions needed to implement kdtree
  integer(ip)                       :: npoib_perset,nboun_perset
  integer(ip), pointer              :: tag(:),npoin_tmp(:),local2global(:),ltypb_perset(:),lnodb_perset(:,:)
  real(rp),    pointer              :: coord_perset(:,:)
  real(rp)                          :: bobox_loc(3,2)
  real(rp),    pointer              :: fabox_loc(:,:,:) => null()
  real(rp),    pointer              :: sabox_loc(:,:,:) => null()
  integer(ip), pointer              :: blink_loc(:)     => null()
  integer(ip), pointer              :: stru2_loc(:)     => null()
  real(rp),    pointer              :: ldist_loc(:)     => null()
  type(netyp), pointer              :: lnele_loc(:)     => null()
  !definitions needed to implement dpopar
  real(rp)                          :: chkdi, vec_proje(ndime), proje(ndime), xcoord(ndime)
  real(rp) :: gap

  allocate(tag(npoin))
  allocate(npoin_tmp(npoin))
  allocate(coord_perset(ndime,npoin))
  allocate(local2global(npoin))
  allocate(ltypb_perset(nboun))
  allocate(lnodb_perset(mnodb,nboun))
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !initiated01: &
  !if( .not.initiated ) then
    !
    if(inotmaster) then
    !mark those nodes (global numeration) which are located in the boundary
    tag(1:npoin) = 0_ip
    do iboun = 1,nboun
      pblty=ltypb(iboun)
      pnodb=nnode(pblty)
      do inodb = 1,pnodb
        tag(lnodb(inodb,iboun)) = 1_ip
      end do
    end do

    !store the coordinates of the boundary nodes in coord_perset
    !npoib_perset is the number of nodes in the boundary
    npoib_perset = 0_ip
    coord_perset = 0.0_rp
    do ipoin = 1,npoin
      if (tag(ipoin) == 1_ip) then
        npoib_perset = npoib_perset + 1_ip
        local2global(npoib_perset) = ipoin !link boundary nodel local and global numbering
        coord_perset(1:ndime,npoib_perset) = coord(1:ndime,ipoin) + displ(1:ndime,ipoin,1)
      end if
    end do

    !npoin_tmp is a local vector which has the global numeration
    npoin_tmp(1:npoin) = 0_ip
    do ipoin = 1,npoib_perset
      npoin_tmp(local2global(ipoin)) = ipoin
    end do

    nboun_perset = 0 !< mistake
    do iboun = 1,nboun
      pblty=ltypb(iboun)
      pnodb=nnode(pblty)
      nboun_perset = nboun_perset + 1_ip
      ltypb_perset(iboun) = ltypb(iboun)
      do inodb = 1,pnodb !< mistake
        lnodb_perset(inodb,iboun) = npoin_tmp(lnodb(inodb,iboun))
      end do
    end do

    call kdtree(&
               1_ip, mnodb, npoib_perset, nboun_perset,&
               coord_perset, lnodb_perset,ltypb_perset,&
               fabox_loc,bobox_loc,sabox_loc,blink_loc,&
               stru2_loc,ldist_loc,lnele_loc)

    do ipoin = 1, n_send
      gap   = huge(1.0_rp)
      chkdi = huge(1.0_rp)  !< mistake
      xcoord(1:ndime) = coords_j(1:ndime,ipoin)
      proje(1:ndime)  = 0.0_rp
      vec_proje(1:ndime)  = 0.0_rp
      call dpopar(&
                 1_ip, xcoord(1:ndime), &
                 npoib_perset, mnodb, nboun_perset,&
                 chkdi, &
                 ltypb_perset, lnodb_perset, coord_perset,&
                 gap, vec_proje, proje, iboun)
      !
      if( present(Dn) ) Dn(        ipoin ) = gap
      if( present(Pn) ) Pn(1:ndime,ipoin ) = proje(1:ndime)
      !
    end do

    !deallocate kdtree
    call kdtree(&
               2_ip,mnodb,npoib_perset,nboun_perset,&
               coord_perset,lnodb_perset,ltypb_perset,&
               fabox_loc,bobox_loc,sabox_loc,blink_loc,&
               stru2_loc,ldist_loc,lnele_loc)
    end if !inotmaster
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  !endif initiated01
  !
  deallocate(tag)
  deallocate(coord_perset)
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_kdtree_01(n_recv, coords_j)
  use mod_kdtree,        only: kdtree, dpopar
  use def_domain,        only: mnodb
  use def_master,        only: netyp, title
  use def_domain,        only: npoib,nboun,lnodb,ltypb
  use def_master,        only: displ
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip), intent(in)            :: n_recv
  real(rp),    intent(inout)         :: coords_j(:)

  integer(ip)                       :: iboun, inodb, ipoin, idime
  integer(ip)                       :: pblty, pnodb

  !definitions needed to implement kdtree
  integer(ip)                       :: npoib_perset,nboun_perset
  integer(ip), pointer              :: tag(:),npoin_tmp(:),local2global(:),ltypb_perset(:),lnodb_perset(:,:)
  real(rp),    pointer              :: coord_perset(:,:)
  real(rp)                          :: bobox_loc(3,2)
  real(rp),    pointer              :: fabox_loc(:,:,:) => null()
  real(rp),    pointer              :: sabox_loc(:,:,:) => null()
  integer(ip), pointer              :: blink_loc(:)     => null()
  integer(ip), pointer              :: stru2_loc(:)     => null()
  real(rp),    pointer              :: ldist_loc(:)     => null()
  type(netyp), pointer              :: lnele_loc(:)     => null()
  !definitions needed to implement dpopar
  real(rp)                          :: chkdi, vec_proje(ndime), proje(ndime), xcoord(ndime)
  real(rp) :: gap

  allocate(tag(npoin))
  allocate(npoin_tmp(npoin))
  allocate(coord_perset(ndime,npoin))
  allocate(local2global(npoin))
  allocate(ltypb_perset(nboun))
  allocate(lnodb_perset(mnodb,nboun))
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !initiated01: &
  !if( .not.initiated ) then
    !
    if(inotmaster) then
    !mark those nodes (global numeration) which are located in the boundary
    tag(1:npoin) = 0_ip
    do iboun = 1,nboun
      pblty=ltypb(iboun)
      pnodb=nnode(pblty)
      do inodb = 1,pnodb
        tag(lnodb(inodb,iboun)) = 1_ip
      end do
    end do

    !store the coordinates of the boundary nodes in coord_perset
    !npoib_perset is the number of nodes in the boundary
    npoib_perset = 0_ip
    coord_perset = 0_rp
    do ipoin = 1,npoin
      if (tag(ipoin) == 1_ip) then
        npoib_perset = npoib_perset + 1_ip
        local2global(npoib_perset) = ipoin !link boundary nodel local and global numbering
        coord_perset(1:ndime,npoib_perset) = coord(1:ndime,ipoin) + displ(1:ndime,ipoin,1)
      end if
    end do

    !npoin_tmp is a local vector which has the global numeration
    npoin_tmp(1:npoin) = 0_ip
    do ipoin = 1,npoib_perset
      npoin_tmp(local2global(ipoin)) = ipoin
    end do

    nboun_perset = 0 !< mistake
    do iboun = 1,nboun
      pblty=ltypb(iboun)
      pnodb=nnode(pblty)
      nboun_perset = nboun_perset + 1_ip
      ltypb_perset(iboun) = ltypb(iboun)
      do inodb = 1,pnodb !< mistake
        lnodb_perset(inodb,iboun) = npoin_tmp(lnodb(inodb,iboun))
      end do
    end do

    call kdtree(&
               1_ip, mnodb, npoib_perset, nboun_perset,&
               coord_perset, lnodb_perset,ltypb_perset,&
               fabox_loc,bobox_loc,sabox_loc,blink_loc,&
               stru2_loc,ldist_loc,lnele_loc)

    do ipoin = 1, n_recv
      gap   = huge(1.0_rp)
      chkdi = huge(1.0_rp)  !< mistake
      proje(1:ndime)  = 0.0_rp
      vec_proje(1:ndime)  = 0.0_rp

      do idime = 1,ndime
        xcoord(idime) = coords_j( ndime*(ipoin-1)+idime )
      enddo

      call dpopar(&
                 1_ip, xcoord(1:ndime), &
                 npoib_perset, mnodb, nboun_perset,&
                 chkdi, &
                 ltypb_perset, lnodb_perset, coord_perset,&
                 gap, vec_proje, proje, iboun)
      !
      do idime = 1,ndime
        coords_j( ndime*(ipoin-1)+idime ) = proje(idime)
      enddo
      !
    end do

    !deallocate kdtree
    call kdtree(&
               2_ip,mnodb,npoib_perset,nboun_perset,&
               coord_perset,lnodb_perset,ltypb_perset,&
               fabox_loc,bobox_loc,sabox_loc,blink_loc,&
               stru2_loc,ldist_loc,lnele_loc)
    end if !inotmaster
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  !endif initiated01
  !
  deallocate(tag)
  deallocate(coord_perset)
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_statistics( n_ij, n_ji )
  use mod_communications, only : PAR_SUM, PAR_MAX, PAR_MIN, PAR_GATHER
  use mod_parall,         only : PAR_COMM_MY_CODE, PAR_CODE_SIZE
  use def_master,         only : title, INOTSLAVE

  implicit none
  integer(ip), intent(in   )   :: n_ij ! PLEPP%n_ij   = n_send
  integer(ip), intent(in   )   :: n_ji ! PLEPP%n_ji   = n_recv
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(4)  :: MPI_RANK, MPI_SIZE, istat4
  integer(ip) :: Nsend, Nrecv, Ncoupled, Npts, Ntrscts
  integer(ip) :: Msend, Mrecv, Mtrscts
  integer(ip) :: ii
  !
  real(rp)    :: Asend, Arecv, Apts, Atrscts
  !
  real(rp)    :: toSend
  real(rp), pointer :: toRecv(:) => null()
  !
  character(len=*), parameter :: FMT2 = '(A, 2I5, 4F12.2)'
  character(100)     :: filename
  integer, parameter :: filehandle = 4321
  !
  !-----------------------------------------------------------------------||---!
  !
if(STATISTICS) then
  MPI_RANK = 1
  MPI_SIZE = 1
#ifndef MPI_OFF
  call MPI_Comm_size(PAR_COMM_MY_CODE, MPI_SIZE, istat4)
  call MPI_Comm_rank(PAR_COMM_MY_CODE, MPI_RANK, istat4)
#endif

  Npts = npoin
  call PAR_SUM(Npts,'IN MY CODE')

  Ncoupled = 0
  if(n_ij > 0) Ncoupled=1 ! coupled subdomain
  call PAR_SUM(Ncoupled,'IN MY CODE')

  Msend = n_ij
  Mrecv = n_ji
  call PAR_MAX(Msend,'IN MY CODE')
  call PAR_MAX(Mrecv,'IN MY CODE')

  Nsend = n_ij
  Nrecv = n_ji
  call PAR_SUM(Nsend,'IN MY CODE')
  call PAR_SUM(Nrecv,'IN MY CODE')

  Ntrscts = 0
#ifdef COMMDOM
  call commdom_locator_get_get_n_intersects( Ntrscts )
#endif
  call PAR_SUM(Ntrscts,'IN MY CODE')
  call PAR_MAX(Mtrscts,'IN MY CODE')

  if(ISEQUEN) then
    Apts    = Npts
    Asend   = Nsend
    Arecv   = Nrecv
    Atrscts = Ntrscts
  else
    Apts    = Npts    / (MPI_SIZE-1)
    Asend   = Nsend   / Ncoupled
    Arecv   = Nrecv   / Ncoupled
    Atrscts = Ntrscts / Ncoupled
  endif

 !if(INOTSLAVE) then
 !  print *, "[commdom_dynamic_statistics] '", trim(title),"' ",  MPI_SIZE, Ncoupled, Asend, Arecv, Atrscts, Apts
 !endif

  if(INOTSLAVE.or.ISEQUEN) then
    if(.not.associated(toRecv)) allocate( toRecv(MPI_SIZE) )
    toRecv   = -1
  else
    if(.not.associated(toRecv)) allocate( toRecv(0) )
  endif

  toSend   = n_ij
  call PAR_GATHER(toSend, toRecv,'IN MY CODE')

  if(INOTSLAVE.or.ISEQUEN) then
    write(*,FMT2) " [commdom_dynamic_statistics] '"//trim(title)//"' ", MPI_SIZE, Ncoupled, Atrscts, Asend, Arecv, Apts
    write(filename,'(a,"_",i6.6,".sta")') trim(title), MPI_SIZE
    open(filehandle, file=trim(filename), STATUS='REPLACE')
    write(filehandle,*) "##MPI_SIZE, Ncoupled, Atrscts, Asend, Arecv, Apts"
   !write(filehandle,*) MPI_SIZE, Ncoupled, Atrscts, Asend, Arecv, Apts
    write(filehandle,'(A, 3I6.6, 3F12.2)') " ##", MPI_SIZE, Ncoupled, Atrscts, minval(toRecv,mask=toRecv>0), Asend, maxval(toRecv) !< 2017JAN09
    do ii=1,MPI_SIZE
      write(filehandle,*) ii, int( toRecv(ii) )
    enddo
    close(filehandle)
  endif

!  deallocate( toRecv )
endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| TO_KILL |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_aitken_relaxation_xxx()
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  ! i: Stage
  ! v: fixed point iteration
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip), parameter :: v = 2
  real(rp), pointer      :: PHI(:,:) !   PHI  = <PHI(v+2), PHI(v+1), PHI(v+0)>
  real(rp), pointer      ::   R(:,:) ! R(v+1) = PHI(v+1) - PHI(v+0)
  !
  real(rp)               :: modR
  real(rp)               :: tolerance
  logical(ip)            :: convergence
  integer(ip)            :: n_inter
  !
  real(rp)               :: omega(3)
  real(rp)               :: omega_0   = -1, omega_n = -1
  real(rp)               :: omega_max = 0.8 !< philipp, tobias, 2014
  !
  R(:,v+1) = PHI(:,v+1) - PHI(:,v+0) ! R(v+1) = PHI(v+1) - PHI(v  )
  R(:,v+0) = PHI(:,v+0) - PHI(:,v-1) ! R(v  ) = PHI(v  ) - PHI(v-1)
  !
  modR        = dot_product( R(:,v+1), R(:,v+1) )
  convergence = sqrt( modR/n_inter ) < tolerance
  !
  ! omega(v+1) = -omega(v) * R(v) * [R(v+1) - R(v)] / [R(v+1) - R(v)]**2
  !
  omega(v+1) = &
               dot_product( R(:,v+0)           , R(:,v+1) - R(:,v+0) ) / &
               dot_product( R(:,v+1) - R(:,v+0), R(:,v+1) - R(:,v+0) )
  !
  omega(v+1) = -omega(v+1) * omega(v)
  !
  !< ulrich wolfgang 2007
  omega_0 = max( omega_n, omega_max)
  !
  !< Dregroote, Souto, 2010
  omega_0 =  min( abs(omega_n), omega_max)
  omega_0 = sign(      omega_0, omega_n  )
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  function commdom_dynamic_init_coupling_00(PLEPP, total, debug) result(ok)
  use mod_communications, only: PAR_SUM
  use def_master,         only: title
  implicit none
  type(COMMDOM_PLEPP_COUPLING), intent(inout) :: PLEPP
  integer(ip),        optional, intent(  out) :: total
  logical(ip),        optional, intent(   in) :: debug
  !
  logical(ip)   :: ok
  integer(ip)   :: n_recv, n_send
  real(rp)      :: prop_sum
  !
  ok       = .false.
  prop_sum = 0.0
  !
  if( (PLEPP%commij /= -1).or.(PLEPP%commji /= -1) ) then
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    n_recv = 0
    n_send = 0
    call commdom_locator_get_n_dist_points( n_send )
    call commdom_locator_get_n_interior(    n_recv )
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    prop_sum = n_send
    call PAR_SUM( prop_sum, 'IN MY CODE')
    !
    if( present(total) ) total = prop_sum
    if( prop_sum /= 0  )    ok = .true.
    if( debug ) print *, "[commdom_dynamic_init_coupling_00] '", trim(title), "' total, partial", prop_sum, n_send
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  endif
  end function
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_dynamic_init_coupling( id_fixbo_j )
  use mod_commdom_plepp,  only: PLEPP_CPLNG
  implicit none
  integer(ip), optional, intent(in) :: id_fixbo_j
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
!    call commdom_dynamic_create(         PLEPP_CPLNG, id_fixbo_j )
    if( commdom_dynamic_init_coupling_00(PLEPP_CPLNG, debug=.true.) ) call runend("EXIT!!")
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

#endif
end module mod_commdom_dynamic
!==============================================================================!
  !-----------------------------------------------------| OUTER_ITERATIONS |---!
  !
  !   |_Alya
  !     |_call Turnon()
  !     |_call Iniunk()
  !     |_time: do while
  !       |_call Timste()
  !       |_reset: do
  !         |_call Begste()                                    coupling_driver_iteration(1:max_block_cou)  = 0
  !           |_block: do while
  !             |_coupling: do while
  !               |_call Begzon()                              coupling_driver_iteration( iblok ) += 1
  !               |_modules: do while
  !                 |_call Doiter()
  !                 |_call Concou()
  !               |_call Endzon()                              call COU_CHECK_CONVERGENCE
  !                                                            call cou_cvgunk
  !             |_call Conblk()                                coupling_driver_iteration( iblok )  = 0
  !       |_call Endste()
  !
  !-----------------------------------------------------------------------||---!
!==============================================================================!
