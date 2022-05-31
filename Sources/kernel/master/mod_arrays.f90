!-----------------------------------------------------------------------
!> @addtogroup Master
!> @{
!> @file    mod_arrays.f90
!> @author  houzeaux
!> @date    2019-11-13
!> @brief   Manipulate arrays
!> @details Manipulate arrays:
!>          1. Register
!>          2. Allocate
!>          3. Deallocate
!>          4. Redistribute
!>          5. Write in restart file
!>          6. Read in restart file
!-----------------------------------------------------------------------

module mod_arrays

  use def_kintyp
  use def_master
  use def_kintyp_mesh_basic, only : mesh_type_basic
  use def_kermod,            only : witness_mesh
  use def_domain,            only : nelem,nboun,npoin
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_size
  use mod_redistribution,    only : redistribution_array
  use mod_postpr,            only : postpr_read_restart
  use mod_postpr,            only : postpr_write_restart
  use mod_postpr,            only : postpr_postprocess
  use mod_communications,    only : PAR_MAX
  use mod_messages,          only : messages_live
  use mod_AMR_interpolate,   only : AMR_interpolate

  implicit none

  integer(8)      :: memor_loc(2)
  integer(ip)     :: enti_posit_loc
  integer(ip)     :: id_loc          ! Mesh ID
  integer(ip)     :: nenti_loc
  character(200)  :: variable_name_loc
  character(200)  :: subroutine_name_loc
  integer(ip)     :: modul_loc
  character(5)    :: wopos_loc(5)
  integer(ip)     :: kfl_reawr_loc
  integer(ip)     :: minmem
  
  character(5)    :: variable_name_5
  integer(8)      :: memor_arrays(2)
  integer(ip)     :: array_total_number(mmodu)
  private

  interface arrays
     module procedure &
          &           arrays_RP_1,arrays_RP_2,arrays_RP_3,arrays_R3P_1,&
          &           arrays_IP_1
  end interface arrays

  interface arrays_used
     module procedure &
          &           arrays_used_RP_1,arrays_used_RP_2,&
          &           arrays_used_RP_3,arrays_used_R3P_1
  end interface arrays_used

  interface arrays_register
     module procedure &
          &           arrays_register_0,&
          &           arrays_register_RP_1,arrays_register_RP_2,&
          &           arrays_register_RP_3,arrays_register_R3P_1,&
          &           arrays_register_IP_1,&
          &           arrays_register_IP_2
  end interface arrays_register

  public :: arrays
  public :: arrays_initialization
  public :: arrays_register
  public :: arrays_number
  public :: arrays_allocated
  public :: arrays_tag             ! wopos of a module variable
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Variable number
  !> @details Find the number of a variable
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function arrays_number(wname,MODULE_NUMBER)

    character(len=*),              intent(in)    :: wname
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    integer(ip)                                  :: ivari,imodu
    
    if( present(MODULE_NUMBER) ) then
       imodu = MODULE_NUMBER
    else
       imodu = modul
    end if

    arrays_number = 0_ip
    do ivari = 1,nvarp
       if( momod(imodu) % postp(1) % wopos (1,ivari) == trim(wname) ) then
          arrays_number = ivari
          return
       end if
    end do
       
  end function arrays_number
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Initialization
  !> @details Initialization of the module
  !> 
  !-----------------------------------------------------------------------
  
  subroutine arrays_initialization()

    memor_arrays       = 0_8
    array_total_number = 0
    minmem             = 0
    
  end subroutine arrays_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Register arrays
  !> @details Register an array:
  !>
  !>          ARRAY_NUMBER ...... Number of the array, must be unique
  !>          WOPOS(1) .......... Name 
  !>          WOPOS(2) .......... Dimension=                  SCALA/VECTO/MATR/R3PVE
  !>          WOPOS(3) .......... Entity=                     NPOIN/NELEM/NBOUN
  !>          WOPOS(4) .......... Primary or secondary array= PRIMA/SECON
  !>          WOPOS(5) .......... Kind of the array=          REA/INTEG
  !>
  !>          ENTITY_POSITION ... Where the entity dimension is located
  !>          TIME_POSITION ..... For transient variables, where time is
  !>
  !>          Primary variables are those required for a restart, a repartitioning
  !>          or a remeshing. Secondary are basically used for postprocess
  !>
  !-----------------------------------------------------------------------
  
  subroutine arrays_register_0(ARRAY_NUMBER,wopos,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER)

    integer(ip),                   intent(in)    :: ARRAY_NUMBER
    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER)

  end subroutine arrays_register_0
  
  subroutine arrays_register_IP_1(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER)

    integer(ip),                   intent(in)    :: ARRAY_NUMBER
    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         pointer,  intent(inout) :: xx(:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER

    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,VAR_INT=.true.)
    
  end subroutine arrays_register_IP_1
  
  subroutine arrays_register_IP_2(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER)

    integer(ip),                   intent(in)    :: ARRAY_NUMBER
    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         pointer,  intent(inout) :: xx(:,:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER

    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 2 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if
    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,VAR_INT=.true.)
    !
    ! Guess component position
    !
    if( .not. present(COMPONENT_POSITION) ) then
       if( wopos(2) == 'SCALA' ) then
          if( .not. present(TIME_POSITION) ) then
             if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 2 ) then
                momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 1
             else if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 1 ) then
                momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 2
             end if
          end if
       end if
    end if
    
  end subroutine arrays_register_IP_2
  
  subroutine arrays_register_RP_1(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER)

    integer(ip),                   intent(in)    :: ARRAY_NUMBER
    character(len=*),              intent(in)    :: wopos(:)
    real(rp),            pointer,  intent(inout) :: xx(:)
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER

    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER)
    
  end subroutine arrays_register_RP_1
  
  subroutine arrays_register_RP_2(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST)

    integer(ip),                             intent(in)    :: ARRAY_NUMBER
    character(len=*),                        intent(in)    :: wopos(:)
    real(rp),                      pointer,  intent(inout) :: xx(:,:)
    integer(ip),         optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),         optional,           intent(in)    :: TIME_POSITION
    integer(ip),         optional,           intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional,           intent(in)    :: EXCLUDE_RST(:)
    integer(ip),         optional,           intent(in)    :: MODULE_NUMBER
    
    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 2 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if
    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST=EXCLUDE_RST)
    !
    ! Guess component position
    !
    if( .not. present(COMPONENT_POSITION) ) then
       if( wopos(2) == 'SCALA' ) then
          if( .not. present(TIME_POSITION) ) then
             if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 2 ) then
                momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 1
             else if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 1 ) then
                momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 2
             end if
          end if
       end if
    end if
    
  end subroutine arrays_register_RP_2
  
  subroutine arrays_register_RP_3(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST)

    integer(ip),                           intent(in)    :: ARRAY_NUMBER
    character(len=*),                      intent(in)    :: wopos(:)
    real(rp),                    pointer,  intent(inout) :: xx(:,:,:)
    integer(ip),       optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),       optional,           intent(in)    :: TIME_POSITION
    integer(ip),       optional,           intent(in)    :: COMPONENT_POSITION
    integer(ip),       optional,           intent(in)    :: MODULE_NUMBER
    integer(ip),       optional,           intent(in)    :: EXCLUDE_RST(:)

    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 3 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if

    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST=EXCLUDE_RST)
    !
    ! Guess component position
    !
    if( .not. present(COMPONENT_POSITION) ) then
       if( wopos(2) == 'SCALA' ) then
          if( present(ENTITY_POSITION) ) then
             if( momod(modul_loc) % postp(1) % time_posit(ARRAY_NUMBER) == 3 ) then
                if(      ENTITY_POSITION == 2 ) then
                   momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 1
                else if( ENTITY_POSITION == 1 ) then
                   momod(modul_loc) % postp(1) % comp_posit(ARRAY_NUMBER) = 2
                end if
             else
                call runend('ARRAYS_REGISTER: CANNOT GUESS')
             end if
          end if
       end if
    end if

  end subroutine arrays_register_RP_3
  
  subroutine arrays_register_R3P_1(ARRAY_NUMBER,wopos,xx,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST)

    integer(ip),                             intent(in)    :: ARRAY_NUMBER
    character(len=*),                        intent(in)    :: wopos(:)
    type(r3p),                     pointer,  intent(inout) :: xx(:)
    integer(ip),         optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),         optional,           intent(in)    :: MODULE_NUMBER
    integer(ip),         optional,           intent(in)    :: TIME_POSITION
    integer(ip),         optional,           intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional,           intent(in)    :: EXCLUDE_RST(:)

    if( present(ENTITY_POSITION) ) then
       if( ENTITY_POSITION > 1 ) call runend('ARRAYS_REGISTER: WRONG POSITION')
    end if
    nullify(xx)
    call arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,EXCLUDE_RST=EXCLUDE_RST)
    
  end subroutine arrays_register_R3P_1
  
  subroutine arrays_register_go(wopos,ARRAY_NUMBER,ENTITY_POSITION,TIME_POSITION,COMPONENT_POSITION,MODULE_NUMBER,VAR_INT,EXCLUDE_RST)

    character(len=*),              intent(in)    :: wopos(:)
    integer(ip),         optional, intent(in)    :: ARRAY_NUMBER
    integer(ip),         optional, intent(in)    :: ENTITY_POSITION
    integer(ip),         optional, intent(in)    :: TIME_POSITION
    integer(ip),         optional, intent(in)    :: COMPONENT_POSITION
    integer(ip),         optional, intent(in)    :: MODULE_NUMBER
    logical(lg),         optional, intent(in)    :: VAR_INT
    integer(ip),         optional, intent(in)    :: EXCLUDE_RST(:)
    integer(ip)                                  :: ii,itime,ivari
    logical(lg)                                  :: if_var_int
    
    if( present(MODULE_NUMBER) ) then
       modul_loc = MODULE_NUMBER
    else
       modul_loc = modul
    end if
    if( present(ARRAY_NUMBER) ) then
       ivari = ARRAY_NUMBER
    else
       array_total_number(modul_loc) = array_total_number(modul_loc) + 1
       ivari = array_total_number(modul_loc)
    end if
    if( present(VAR_INT) ) then
       if_var_int = VAR_INT
    else
       if_var_int = .false.       
    end if
    if( present(EXCLUDE_RST) ) then
       do ii = 1,size(EXCLUDE_RST)
          itime = EXCLUDE_RST(ii)
          momod(modul_loc) % postp(1) % rst_time(itime,ivari) = .false. 
       end do
    end if
    
    if( momod(modul_loc) % postp(1) % array_registered(ivari) == 1 ) then
       call runend('ARRAY_REGISTER: ARRAY '//trim(wopos(1))//' HAS ALREADY BEEN REGISTERED')
    else
       momod(modul_loc) % postp(1) % array_registered(ivari) = 1
       momod(modul_loc) % postp(1) % wopos (1,ivari)         = trim(wopos(1))
       momod(modul_loc) % postp(1) % wopos (2,ivari)         = trim(wopos(2))
       momod(modul_loc) % postp(1) % wopos (3,ivari)         = trim(wopos(3))
       momod(modul_loc) % postp(1) % wopos (4,ivari)         = trim(wopos(4))
       if( size(wopos) < 5 ) then
          if( if_var_int ) then
             momod(modul_loc) % postp(1) % wopos (5,ivari)      = 'INT'
          else
             momod(modul_loc) % postp(1) % wopos (5,ivari)      = 'REAL'
          end if
       else
          momod(modul_loc) % postp(1) % wopos (5,ivari)      = trim(wopos(5))
       end if
       
       if( present(ENTITY_POSITION) ) then
          momod(modul_loc) % postp(1) % enti_posit (ivari) = ENTITY_POSITION
       else       
          momod(modul_loc) % postp(1) % enti_posit (ivari) = 1
       end if
       
       if( present(TIME_POSITION) ) then
          momod(modul_loc) % postp(1) % time_posit (ivari) = TIME_POSITION
       else       
          momod(modul_loc) % postp(1) % time_posit (ivari) = 0
       end if
       
       if( present(COMPONENT_POSITION) ) then
          momod(modul_loc) % postp(1) % comp_posit (ivari) = COMPONENT_POSITION
       else       
          momod(modul_loc) % postp(1) % comp_posit (ivari) = 0
       end if
    end if
 
  end subroutine arrays_register_go

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Do what we have to do with arrays
  !> @details Do what we have to do with arrays
  !> 
  !-----------------------------------------------------------------------

  subroutine arrays_RP_1(ivari,wtask,xx,ndim1,ENTITY_POSITION,TIME_POSITION,TAG1,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari
    character(len=*),                           intent(in)    :: wtask
    real(rp),                         pointer,  intent(inout) :: xx(:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    integer(ip),            optional,           intent(in)    :: TAG1
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID
    
    call arrays_options_start(ivari)
    
    select case ( trim(wtask) )
       
    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) ) then
          if( enti_posit_loc == 1 .and. ndim1 < nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc /= 1 )                          call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( trim(wopos_loc(4)) /= 'PRIMA' )           call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1))
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if
       
    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0
       
    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_write_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
         
    case ( 'READ RESTART' )
       !
       ! Write restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_read_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
      end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !  
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
      end if
       
    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))
       
    end select

    call arrays_options_end()

  end subroutine arrays_RP_1
  
  subroutine arrays_RP_2(ivari,wtask,xx,ndim1,ndim2,ENTITY_POSITION,TIME_POSITION,TAG1,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari
    character(len=*),                           intent(in)    :: wtask
    real(rp),                         pointer,  intent(inout) :: xx(:,:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ndim2
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    integer(ip),            optional,           intent(in)    :: TAG1
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID
    
    call arrays_options_start(ivari)

    select case ( trim(wtask) )

    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) .and. present(ndim2) ) then
          if( enti_posit_loc == 1 .and. ndim1 < nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc == 2 .and. ndim2 < nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( trim(wopos_loc(4)) /= 'PRIMA' )           call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1),max(minmem,ndim2))
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
          if(      postp(1) % time_posit(ivari) == 1 ) then
             postp(1) % time_num(ivari) = ndim1
          else if( postp(1) % time_posit(ivari) == 2 ) then
             postp(1) % time_num(ivari) = ndim2
          end if
          if(      postp(1) % comp_posit(ivari) == 1 ) then
             postp(1) % comp_num(ivari) = ndim1
          else if( postp(1) % comp_posit(ivari) == 2 ) then
             postp(1) % comp_num(ivari) = ndim2
          end if
          if( postp(1) % wopos(2,ivari) == 'VECTO' ) then
             if(      postp(1) % comp_posit(ivari) /= 1 .and. postp(1) % time_posit(ivari) /= 1 ) then
                postp(1) % dime_num(ivari) = ndim1
             else 
                postp(1) % dime_num(ivari) = ndim2
             end if
          end if
          
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if
       
    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0
       
    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)

    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_write_restart(xx,ivari,ittim,cutim,posit=enti_posit_loc,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'READ RESTART' )
       !
       ! Read restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_read_restart(xx,ivari,ittim,cutim,posit=enti_posit_loc,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
      end if

    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))
       
    end select

    call arrays_options_end()

  end subroutine arrays_RP_2
  
  subroutine arrays_RP_3(ivari,wtask,xx,ndim1,ndim2,ndim3,ENTITY_POSITION,TIME_POSITION,TAG1,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari    
    character(len=*),                           intent(in)    :: wtask
    real(rp),                         pointer,  intent(inout) :: xx(:,:,:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ndim2
    integer(ip),            optional,           intent(in)    :: ndim3
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    integer(ip),            optional,           intent(in)    :: TAG1
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID

    call arrays_options_start(ivari)

    select case ( trim(wtask) )

    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) .and. present(ndim2) .and. present(ndim3) ) then
          if( enti_posit_loc == 1 .and. ndim1 /= nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc == 2 .and. ndim2 /= nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc == 3 .and. ndim3 /= nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( trim(wopos_loc(4)) /= 'PRIMA' )           call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1),max(minmem,ndim2),max(minmem,ndim3))
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
          if(      postp(1) % time_posit(ivari) == 1 ) then
             postp(1) % time_num(ivari) = ndim1
          else if( postp(1) % time_posit(ivari) == 2 ) then
             postp(1) % time_num(ivari) = ndim2
          else if( postp(1) % time_posit(ivari) == 3 ) then
             postp(1) % time_num(ivari) = ndim3
          end if
          if(      postp(1) % comp_posit(ivari) == 1 ) then
             postp(1) % comp_num(ivari) = ndim1
          else if( postp(1) % comp_posit(ivari) == 2 ) then
             postp(1) % comp_num(ivari) = ndim2
          else if( postp(1) % comp_posit(ivari) == 3 ) then
             postp(1) % comp_num(ivari) = ndim3
          end if
          if( postp(1) % wopos(2,ivari) == 'VECTO' ) then
             if(      postp(1) % comp_posit(ivari) /= 1 .and. postp(1) % time_posit(ivari) /= 1 ) then
                postp(1) % dime_num(ivari) = ndim1
             else if( postp(1) % comp_posit(ivari) /= 2 .and. postp(1) % time_posit(ivari) /= 2 ) then
                postp(1) % dime_num(ivari) = ndim2
             else
                postp(1) % dime_num(ivari) = ndim3
             end if
          end if          
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if

    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0

    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)

    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)

    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')          
          call postpr_write_restart(xx,ivari,ittim,cutim,TAG1=TAG1)
          call messages_live('VARIABLE','END SECTION')
       end if

    case ( 'READ RESTART' )
       !
       ! Read restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')          
          call postpr_read_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
      end if

    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))

    end select

    call arrays_options_end()

  end subroutine arrays_RP_3
  
  subroutine arrays_IP_1(ivari,wtask,xx,ndim1,ENTITY_POSITION,TIME_POSITION,TAG1,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari
    character(len=*),                           intent(in)    :: wtask
    integer(ip),                      pointer,  intent(inout) :: xx(:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    integer(ip),            optional,           intent(in)    :: TAG1
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID

    call arrays_options_start(ivari)
    
    select case ( trim(wtask) )
       
    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) ) then
          if( enti_posit_loc == 1 .and. ndim1 < nenti_loc ) call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( enti_posit_loc /= 1 )                          call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' NOT REGISTERED CORRECTLY')
          if( trim(wopos_loc(4)) /= 'PRIMA' )           call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1))
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if
       
    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0
       
    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)

    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_write_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
         
    case ( 'READ RESTART' )
       !
       ! Write restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_read_restart(xx,ivari,ittim,cutim,TAG1=TAG1)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !  
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
      end if

    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))
       
    end select

    call arrays_options_end()

  end subroutine arrays_IP_1
    
  subroutine arrays_R3P_1(ivari,wtask,xx,ndim1,ENTITY_POSITION,TIME_POSITION,FORCE,MESH,MESH_ID)

    integer(ip),                                intent(in)    :: ivari
    character(len=*),                           intent(in)    :: wtask
    type(r3p),                        pointer,  intent(inout) :: xx(:)
    integer(ip),            optional,           intent(in)    :: ndim1
    integer(ip),            optional,           intent(in)    :: ENTITY_POSITION
    integer(ip),            optional,           intent(in)    :: TIME_POSITION
    logical(lg),            optional,           intent(in)    :: FORCE
    type(mesh_type_basic),  optional,           intent(in)    :: MESH
    integer(ip),            optional,           intent(in)    :: MESH_ID

    call arrays_options_start(ivari)
    
    select case ( trim(wtask) )
       
    case ( 'ALLOCATE' )
       !
       ! Allocate
       !
       if( present(ndim1) ) then
          if( trim(wopos_loc(4)) /= 'PRIMA' )           call runend('ARRAYS: VARIABLE '//trim(wopos_loc(1))//' SHOULD BE REGISTERED AS PRIMARY')
          call memory_alloca(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx,max(minmem,ndim1))
          postp(1) % array_used(ivari) = 1
          if( memory_size(xx) > 0 ) postp(1) % array_allocated(ivari) = 1 
       else
          call runend('ARRAYS: DIMENSIONS ARE MISSING')
       end if
       
    case ( 'DEALLOCATE' )
       !
       ! deallocate
       !
       call memory_deallo(memor_loc,trim(variable_name_loc),trim(subroutine_name_loc),xx)
       postp(1) % array_used(ivari) = 0
       postp(1) % array_allocated(ivari) = 0
       
    case ( 'REDISTRIBUTE' )
       !
       ! Redistribute
       !
       call redistribution_array(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'INTERPOLATE' )
       ! 
       ! Interpolate
       !
       call runend('NOT CODED')
       !call AMR_interpolate(xx,wopos_loc(3),POSIT=enti_posit_loc,MEMOR=memor_loc,VARIABLE_NAME=variable_name_loc)
       
    case ( 'WRITE RESTART' )
       !
       ! Write restart
       !
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')         
          call postpr_write_restart(xx,ivari,ittim,cutim)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'READ RESTART' )
       !
       ! Write restart
       !  
       if( postp(1) % array_used(ivari) == 1 ) then
          call messages_live('VARIABLE: '//wopos_loc(1)//' FROM '//namod(modul),'START SECTION')
          call postpr_read_restart(xx,ivari,ittim,cutim)   
          call messages_live('VARIABLE','END SECTION')
       end if
       
    case ( 'POSTPROCESS' )
       !
       ! Postprocess
       !  
       if( postp(1) % array_used(ivari) == 1 .or. optional_argument(.false.,FORCE) ) then
          if( present(MESH) ) then
             call runend('arrays_R3P_1: POSTPROCESS NOT CODED 1')                
             !call postpr_postprocess(xx,ivari,ittim,cutim,MESH=MESH)             
          else if( present(MESH_ID) ) then
             if(     MESH_ID == 0 ) then
                call postpr_postprocess(xx,ivari,ittim,cutim)
             else if( MESH_ID > 0 ) then
                call runend('arrays_R3P_1: POSTPROCESS NOT CODED 2')                
                !call postpr_postprocess(xx,ivari,ittim,cutim,MESH=witness_mesh(MESH_ID) % mesh)
             end if
          else
             call postpr_postprocess(xx,ivari,ittim,cutim)
          end if
       end if
       
    case default
       !
       ! Do not know what to do
       !
       call runend('ARRAYS: DO NOT KNOW WHAT TO DO '//trim(wtask)//' WITH '//trim(wopos_loc(1)))
       
    end select

    call arrays_options_end()

  end subroutine arrays_R3P_1
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Starting operations
  !> @details Starting operations to deal with optional arguments
  !> 
  !-----------------------------------------------------------------------

  subroutine arrays_options_start(ivari)
    
    integer(ip),  intent(in) :: ivari

    if( modul == 0 ) then
       call runend('ARRAYS: YOU ARE TRYING TO REGISTER A VARIABLE OUT OF A MODULE')
    else
       kfl_reawr_loc       = kfl_reawr
       modul_loc           = modul
       variable_name_loc   = momod(modul_loc) % postp(1) % wopos(1,ivari)
       wopos_loc           = momod(modul_loc) % postp(1) % wopos(:,ivari)
       subroutine_name_loc = trim(exmod(modul_loc))//'_ARRAYS'
       memor_loc           = mem_modul(1:2,modul_loc)
       enti_posit_loc      = momod(modul_loc) % postp(1) % enti_posit (ivari)
       
       if(      momod(modul_loc) % postp(1) % wopos(3,ivari) == 'NPOIN' ) then
          nenti_loc = npoin
       else if( momod(modul_loc) % postp(1) % wopos(3,ivari) == 'NELEM' ) then
          nenti_loc = nelem
       else if( momod(modul_loc) % postp(1) % wopos(3,ivari) == 'NBOUN' ) then
          nenti_loc = nboun
       end if
    end if
    
  end subroutine arrays_options_start

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-13
  !> @brief   Ending operations
  !> @details Ending operations to deal with optional arguments
  !> 
  !-----------------------------------------------------------------------

  subroutine arrays_options_end()
        
    mem_modul(1:2,modul_loc) = memor_loc
    kfl_reawr = kfl_reawr_loc
    
  end subroutine arrays_options_end

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-15
  !> @brief   Check if variable is allocated
  !> @details Check if variable is allocated at least in one MPI
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function arrays_allocated(wname,MODULE_NUMBER)
    character(len=*),           intent(in) :: wname
    integer(ip),      optional, intent(in) :: MODULE_NUMBER
    integer(ip)                            :: ivari

    ivari = arrays_number(wname,MODULE_NUMBER)

    if( ivari == 0 ) call runend('MOD_ARRAYS: YOU ARE CHECKING A NON REGISTERED VARIABLE '//trim(wname))
    arrays_allocated = .false.
    if( present(MODULE_NUMBER) ) then
       if( momod(MODULE_NUMBER) % postp(1) % array_allocated(ivari) /= 0 )  arrays_allocated = .true.
    else
       if( postp(1) % array_allocated(ivari) /= 0 ) arrays_allocated = .true.
    end if
    
  end function arrays_allocated

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Return array name
  !> @details ???
  !> 
  !-----------------------------------------------------------------------

  function arrays_tag(wname,MODULE_NUMBER) result(wopos)

    character(len=*),              intent(in) :: wname
    integer(ip),         optional, intent(in) :: MODULE_NUMBER
    character(5)                              :: wopos(5)
    integer(ip)                               :: ivari,imodu
    
    if( present(MODULE_NUMBER) ) then
       imodu = MODULE_NUMBER
    else
       imodu = modul
    end if

    wopos = ''
    do ivari = 1,nvarp
       if( momod(imodu) % postp(1) % wopos (1,ivari) == trim(wname) ) then
          wopos(1:5) = momod(imodu) % postp(1) % wopos (1:5,ivari)
          return
       end if
    end do
       
  end function arrays_tag
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-15
  !> @brief   Check if variable is allocated
  !> @details Check if variable is allocated at least in one MPI
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function arrays_used_RP_1(xx)

    real(rp),   pointer, intent(in) :: xx(:)
    integer(ip)                     :: if_allocated

    if_allocated = memory_size(xx)
    call PAR_MAX(if_allocated)
    if( if_allocated > 0 ) then
       arrays_used_RP_1 = .true.
    else
       arrays_used_RP_1 = .false.
    end if
    
  end function arrays_used_RP_1
  
  logical(lg) function arrays_used_RP_2(xx)

    real(rp),   pointer, intent(in) :: xx(:,:)
    integer(ip)                     :: if_allocated

    if_allocated = memory_size(xx)
    call PAR_MAX(if_allocated)
    if( if_allocated > 0 ) then
       arrays_used_RP_2 = .true.
    else
       arrays_used_RP_2 = .false.
    end if
    
  end function arrays_used_RP_2
  
  logical(lg) function arrays_used_RP_3(xx)

    real(rp),   pointer, intent(in) :: xx(:,:,:)
    integer(ip)                     :: if_allocated

    if_allocated = memory_size(xx)
    call PAR_MAX(if_allocated)
    if( if_allocated > 0 ) then
       arrays_used_RP_3 = .true.
    else
       arrays_used_RP_3 = .false.
    end if
    
  end function arrays_used_RP_3
  
  logical(lg) function arrays_used_R3P_1(xx)

    type(r3p),  pointer, intent(in) :: xx(:)
    integer(ip)                     :: if_allocated

    if_allocated = memory_size(xx)
    call PAR_MAX(if_allocated)
    if( if_allocated > 0 ) then
       arrays_used_R3P_1 = .true.
    else
       arrays_used_R3P_1 = .false.
    end if
    
  end function arrays_used_R3P_1
  
end module mod_arrays
!> @}
