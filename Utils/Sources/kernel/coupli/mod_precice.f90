!===================================================================================!
!
!< 2015May25 -> created   
!< 2015Jun02 -> precice_on_field 
!< 
!< 
!
!===================================================================================!
module mod_precice
  use def_kintyp
  !usa def_parall, only : nproc_par, iproc_par
  use def_master, only : inotmaster, imaster, isequen
  use def_master, only : title, dtime   
  use def_domain, only : coord 
  use def_domain, only : npoin, ndime  
  implicit none 
 !integer(ip)   ::  kfl_required ! 1: yes, 0: no

  ! preCICE variables
  integer(ip)            :: &
       kfl_preci,           &  ! preCICE used (=1) or not (=0)
       pre_excha,           &  ! preCICE exchange, i.e. proc on wet surface (=1) or not (=0)
       pre_final               ! =1 if preCICE wants to end the simulation

  integer(ip)            :: &
       pre_idMes,           &  ! preCICE mesh ID
       pre_idF,             &  ! preCICE data ID for forces
       pre_idD,             &  ! preCICE data ID for displacements
       pre_idDD,            &  ! preCICE data ID for displacement deltas
       pre_codno,           &  ! preCICE boundary code
       pre_nvert,           &  ! preCICE number of vertices on wet surface
       pre_ndoub               ! preCICE number of double boundary nodes
  real(rp)               :: &
       pre_tsl                 ! preCICE time step length

  integer(ip), pointer   :: &
       pre_readIDs(:),      & ! preCICE read position IDs
       pre_writeIDs(:),     & ! preCICE write position IDs
       pre_vertIDs(:),      & ! preCICE vertex position IDs
       pre_indre(:),        & ! preCICE read look-up vector for nodes on wet surface
       pre_index(:),        & ! preCICE vertex look-up vector for nodes on wet surface
       pre_indwr(:)           ! preCICE write look-up vector for nodes on wet surface

  real(rp), pointer      :: &
       pre_force(:),        &  ! preCICE forces
       pre_olddi(:)            ! preCICE displ from old fsi iteration
  !
  ! preCICE data structures
  !  
  real(rp),    pointer       :: positions(:)  
  real(rp),    pointer       ::      preD(:) 
  real(rp),    pointer       ::     preDD(:) 
  !
  type     PRECICE_COUPLING
    character(16_ip)           :: name_mesh       = '' 
    integer(ip)                ::   id_mesh       = -1_ip
    !
    integer(ip)                ::    n_props      =  -1_ip 
    integer(ip),      pointer  ::   id_props(:)   => null() 
    integer(ip),      pointer  ::  dim_props(:)   => null()
    integer(ip),      pointer  :: size_props(:)   => null()
    character(16_ip), pointer  :: name_props(:)   => null() 
    !
    real(rp),    pointer       ::      props(:,:)  => null()
    real(rp),    pointer       ::  positions(:)    => null() 
    integer(ip), pointer       :: vertex_ids(:)    => null() 
    !
    integer(ip)                :: idx_vertices_j  = -1_ip      ! 'pre_codno' 
    integer(ip)                ::   n_vertices_j  = -1_ip      ! 'pre_nvert' preCICE number of vertices on wet surface
    integer(ip), pointer       :: idx_coords_j(:) => null()
    !
    real(rp)                   :: tsl             = -huge(1.0) ! 'pre_tsl' preCICE time step length
    !
    logical(ip)                :: wet_processor   = .false.    ! 'pre_excha' is processor not on WetSurface? 
    !
    integer(ip)                :: current_code    = -1_ip   
    integer(ip)                :: code_i          = -1_ip
    integer(ip)                :: code_j          = -1_ip
    integer(ip)                :: module_i        = -1_ip
    integer(ip)                :: module_j        = -1_ip
    !
    logical(ip)                :: debug           = .true.
    !
  end type PRECICE_COUPLING 
  ! 
 !type(PRECICE_COUPLING) PRE_CPLNG 
  ! 
  !
  private
  !
  public :: PRECICE_COUPLING 
  public :: precice_create ! configure  
  public :: precice_initialize
  public :: precice_advance  
  public :: precice_finalize
  public :: precice_sendrecv
  public :: precice_set_time
  !
contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine precice_create( CPLNG ) 
  use def_master, only : current_code
  use def_master, only : ID_TEMPER
  implicit none
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  !
  character(50)   ::  corsc
  integer(ip)     ::  kfl_required ! 1: yes, 0: no
  !
  corsc(:) = ' '
  !-----------------------------------------------------------------------||---!
  !
  ! Create preCICE, called by sld_turnon
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE
    !
    !      <participant name="NASTIN">
    !
    !call precicef_create( trim(title)  , "precice.xml", iproc_par, nproc_par)
    call precicef_action_read_sim_checkp(corsc)
    call precicef_action_required(corsc, kfl_required)
    !
    ! restart data handled by alya directly
    if(kfl_required==1) call precicef_fulfilled_action(corsc)
    !
    !
    CPLNG%n_props        = 1  
    !
    if( .not.associated(CPLNG%id_props)    ) allocate( CPLNG%id_props(CPLNG%n_props)   ) 
    if( .not.associated(CPLNG%dim_props)   ) allocate( CPLNG%dim_props(CPLNG%n_props)  )
    if( .not.associated(CPLNG%size_props)  ) allocate( CPLNG%size_props(CPLNG%n_props) )
    if( .not.associated(CPLNG%name_props)  ) allocate( CPLNG%name_props(CPLNG%n_props) ) 
    !
    if(trim(title) == "DIRICH") then 
      CPLNG%code_i         =  1
      CPLNG%name_mesh      = "Mesh1"
      !
      CPLNG%id_props(1)    = -1_ip  
      CPLNG%dim_props(1)   =  1_ip  
      CPLNG%name_props(1)  = 'TEMPE'
      CPLNG%idx_vertices_j =  1_ip 
      !
      CPLNG%module_i       =  ID_TEMPER 
    endif 
    !
    if(trim(title) == "NEUMAN") then  
      CPLNG%code_j         =  2
      CPLNG%name_mesh      = "Mesh2"
      !
      CPLNG%id_props(1)    = -1_ip
      CPLNG%dim_props(1)   =  1_ip 
      CPLNG%name_props(1)  = 'TFLUX'
      CPLNG%idx_vertices_j =  1_ip
      !
      CPLNG%module_j       =  ID_TEMPER 
    endif 
    !
    CPLNG%current_code     = current_code 
    !
    if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then
      print *, "[precice_init] ", "'", trim(title), ".", trim(CPLNG%name_mesh), "'" 
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
  subroutine precice_initialize( CPLNG )
  use def_domain, only : ltype, lnods, coord, nelem 
  implicit none
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  !
  integer(ip)          :: id 
  !-----------------------------------------------------------------------||---!
  !
  ! Create preCICE, called by sld_iniunk  
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE
    !
    !      <mesh name="WetSurface">
    !
    call precicef_get_mesh_id( trim(CPLNG%name_mesh), CPLNG%id_mesh)
    !
    ! olver-interface dimensions="2">
    !       <data:vector name="Forces"/>
    !
    do id = 1,CPLNG%n_props
      call precicef_get_data_id( trim(CPLNG%name_props(id)), CPLNG%id_mesh, CPLNG%id_props(id) )
    enddo
    ! 
    call precice_on_field(CPLNG, CPLNG%idx_vertices_j) 
    !
    ! pre_nvert   = CPLNG%n_vertices_j, 
    ! pre_index   = CPLNG%idx_coords_j(:)  
    ! pre_vertIDs = CPLNG%vertex_ids(:)     
    !
    if( .not.associated(CPLNG%vertex_ids) ) allocate( CPLNG%vertex_ids(CPLNG%n_vertices_j) ) 
    if( .not.associated(CPLNG%props) ) then  
      do id = 1,CPLNG%n_props
        CPLNG%size_props(id) = CPLNG%dim_props(id) * CPLNG%n_vertices_j 
        allocate( CPLNG%props( CPLNG%size_props(id), id) )
      enddo
    endif  
    ! 
    CPLNG%vertex_ids(:) =  -1_ip
    CPLNG%props(:,:)    = 0.0_rp   
    !
    call precicef_set_vertices(CPLNG%id_mesh, CPLNG%n_vertices_j, CPLNG%positions(1:ndime*CPLNG%n_vertices_j), CPLNG%vertex_ids(1:CPLNG%n_vertices_j) )
    if( associated(CPLNG%positions) ) deallocate( CPLNG%positions )   
    !
    call precicef_initialize( CPLNG%tsl )
    !
    if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then
      print *, "[precice_initialize] ", "'"//trim(title)//"'", CPLNG%tsl  
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
  subroutine precice_advance( CPLNG )
  implicit none
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  !
  character(50)  :: coric, cowic, corsc, cowsc ! constant read/write-iteration/sim-checkpoint
  integer(ip)    :: kfl_required ! 1: yes, 0: no
  !
  cowic(1:50)='                                                  '
  coric(1:50)='                                                  '
  cowsc(1:50)='                                                  '
  corsc(1:50)='                                                  '
  !
  !-----------------------------------------------------------------------||---!
  !
  ! Create preCICE, called by  sld_concou  
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE
  call precicef_action_write_iter_checkp( cowic )
  call precicef_action_required(cowic, kfl_required)
  !
  if(kfl_required == 1) then 
    call precicef_fulfilled_action( cowic )
    kfl_required = -1 
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
!  call precicef_write_bvdata(pre_idDD, pre_nvert, pre_vertIDs, preDD)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call precicef_action_write_sim_checkp( cowsc )
  call precicef_action_required( cowsc, kfl_required )
  if (kfl_required==1) then
    call precicef_fulfilled_action( cowsc )
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call precicef_action_read_iter_checkp( coric )
  call precicef_action_required( coric, kfl_required )
  call precicef_ongoing(kfl_required)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
!
!        pre_tsl = dtime
!        call precicef_advance(pre_tsl)
!
  if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then
    print *, "[precice_advance] ", "'"//trim(title)//"'", "1.0"
  endif
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine precice_convergence( CPLNG )
  use def_master,  only: kfl_gocou
  use def_master,  only: kfl_reset
  implicit none 
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  !
  character(50)  :: coric
  integer(ip)    :: kfl_required ! 1: yes, 0: no
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
        call precicef_action_required(coric, kfl_required)
        !
        if (kfl_required == 1) then 
         !print *, "PRECICE not converged, reset"
          !
          kfl_gocou = 0 !just to be sure
          kfl_reset = 1 !here we use reset
          call precicef_fulfilled_action( coric )
          kfl_required = -1
          !
        else !i.e. converged
         !print *, "PRECICE converged, advancing in time"
          kfl_reset = 0
!          do iinde = 1,pre_nvert
!            do idime = 1,ndime
!              ipoin = pre_index(iinde)
!              pre_olddi((iinde-1)*ndime+idime) = displ(idime,ipoin,1)
!            end do
!          end do
          !
        end if
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
        call precicef_ongoing(kfl_required)
        if (kfl_required==0) then
          print *, "PRECICE asks to shut down"
          pre_final = 1
        end if
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then
    endif 
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine precice_finalize( CPLNG )
  implicit none
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  !
  character(50) :: cowic
  integer(ip)   :: kfl_required ! 1: yes, 0: no
  !
  !-----------------------------------------------------------------------||---!
  !
  ! Create preCICE, called by sld_turnof 
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE
  call precicef_action_write_iter_checkp( cowic )
  call precicef_action_required( cowic, kfl_required )
  !
  if(pre_excha == 1) then
    if( associated(CPLNG%props) ) deallocate( CPLNG%props ) 
  end if
  !
  if(kfl_required == 1) then
    call precicef_fulfilled_action(cowic)
    kfl_required = -1
  endif
  !
  call precicef_finalize()
  !
  if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then
    print *, "[precice_finalize] ", "'"//trim(title)//"'", "1.0"
  endif
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine precice_set_time( CPLNG )
  implicit none 
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  !-----------------------------------------------------------------------||---!
  !
  ! Create preCICE, called by sld_timste  
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE
  dtime = min(dtime, CPLNG%tsl)
  !
  if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then
    print *, "(precice_set_time) ", "'"//trim(title)//"'", "1.0"
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
  subroutine precice_sendrecv( CPLNG, id, to_send, to_recv ) 
  implicit none
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  integer(ip),                  intent(in   )   :: id
  !
  real(rp),  optional,          intent(in   )   :: to_send(npoin)
  real(rp),  optional,          intent(out  )   :: to_recv(npoin)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE
  !
  ! write to preCICE
  if( present(to_send) ) then 
    if(imaster) then 
      CPLNG%props(1:CPLNG%size_props(id), id) = to_send( CPLNG%idx_coords_j(1:CPLNG%n_vertices_j) )  
      call precicef_write_bvdata(id, CPLNG%n_vertices_j, CPLNG%vertex_ids(1:CPLNG%n_vertices_j), CPLNG%props(1:CPLNG%size_props(id), id) )
    endif 
  endif 
  !
  ! read from preCICE
  if( present(to_recv) ) then
    if(imaster) then
      call precicef_read_bvdata(id, CPLNG%n_vertices_j, CPLNG%vertex_ids(1:CPLNG%n_vertices_j), CPLNG%props(id,1:CPLNG%n_vertices_j) )
      to_recv( CPLNG%idx_coords_j(1:CPLNG%n_vertices_j) ) = CPLNG%props(1:CPLNG%size_props(id), id) 
    endif
  endif
  !
  if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then  
    print *, "[precice_sendrecv] ", "'"//trim(title)//"'", "'", trim( CPLNG%name_props(id) ), "'"   
  endif
  !
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine precice_exchange02( CPLNG )
  implicit none
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE 
  ! 
  ! this is not most important call, here everything is done, interpolation, cplscheme, 
  ! communication, output, etc ... whatever you define in the precice-config.xml
  ! the timestep is input and output, you have to tell precice which timestep you want to
  ! use and precice tells you then which timestep you should actually use (for subcycling, 
  ! to sync with other participant etc)
  call precicef_advance( CPLNG%tsl  )
  !
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine precice_checkpoint( CPLNG )
  implicit none
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE 
!  !
!  ! Create
!        call precicef_action_read_sim_checkp(CORSC)
!        call precicef_action_required(CORSC, kfl_required)
!        if (kfl_required==1) then
!          call precicef_fulfilled_action(CORSC)
!        end if
!
!  ! advance 
!        call precicef_action_write_iter_checkp(COWIC)
!        call precicef_action_required(COWIC,i kfl_required)
!        if (kfl_required==1) then 
!          call precicef_fulfilled_action(COWIC)
!          kfl_required = -1
!        end if 
!  !
!        ! no you have to look if you have converged in this timestep
!        ! read iteration checkpoint necessary = no convergence, you should read your back-up   
!        call precicef_action_read_iter_checkp(CORIC) ! get the constant string
!        call precicef_action_required(CORIC, kfl_required)
!        if (kfl_required==1) then !i.e. not yet converged
!          print *, "PRECICE, not converged"
!          kfl_gocou=1 !here we use the coupling loop
!          call precicef_fulfilled_action(CORIC)
!          kfl_required = -1
!  !
!  ! finalize 
!        call precicef_action_write_iter_checkp(COWIC)
!        call precicef_action_required(COWIC, kfl_required)
!        if (kfl_required==1) then !i.e. checkpointing required
!          call precicef_fulfilled_action(COWIC)
!          kfl_required = -1
!        end if
!
!  !
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
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
  subroutine precice_on_field(CPLNG, where_number)
  use def_domain, only:  nfiel, kfl_field, xfiel
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
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
!          print*, coord(:,ipoin), xfiel(ifiel) % a(idime,ipoin)  
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
    if( .not.associated(CPLNG%idx_coords_j) ) allocate( CPLNG%idx_coords_j(CPLNG%n_vertices_j) )
    if( .not.associated(CPLNG%positions)    ) allocate( CPLNG%positions(ndime*CPLNG%n_vertices_j) )
    CPLNG%idx_coords_j = 0.0_rp
    CPLNG%positions    = 0.0_rp
!
!    !
!    if(INOTMASTER) then
!      if(associated(coords_j)) then
!         nullify( coords_j )
!        allocate( coords_j(ndime,CPLNG%n_vertices_j) )
!      endif
!    endif
!    !
!    n_coords_j = CPLNG%n_vertices_j
!
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
         CPLNG%idx_coords_j(i_coords_j) = ipoin
         do idime = 1,ndime 
           dummy = dummy + 1
           CPLNG%positions(dummy) = coord(idime,ipoin)  
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
  if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) print *, "[precice_on_field]", " n_coords_j:", CPLNG%n_vertices_j, "where:", where_number
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine precice_xxx( CPLNG )
  implicit none
  type(PRECICE_COUPLING),       intent(inout)   :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef PRECICE

  if( (INOTMASTER.or.ISEQUEN).and.CPLNG%debug ) then 
    print *, "[precice_xxx]"
  endif 
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
end module mod_precice
