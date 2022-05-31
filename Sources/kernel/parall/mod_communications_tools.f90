!------------------------------------------------------------------------
!> @defgroup Communication_Toolbox
!> Toolbox for MPI communication, bridge to MPI
!> @{
!> @name    Parallelization toolbox
!> @file    mod_communications_tools.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for parallel communications
!> @details ToolBox for parallel communications
!------------------------------------------------------------------------

module mod_communications_tools

  use def_communications
  implicit none
  
#ifdef EXTRAE
  use extrae_module
#endif

  private
  !
  ! Communicator operations
  !
  interface PAR_COMM_RANK_AND_SIZE
     module procedure PAR_COMM_RANK_AND_SIZE_4,PAR_COMM_RANK_AND_SIZE_4W,&
          &           PAR_COMM_RANK_AND_SIZE_41W,&
          &           PAR_COMM_RANK_AND_SIZE_8,PAR_COMM_RANK_AND_SIZE_8W
  end interface PAR_COMM_RANK_AND_SIZE
  !
  ! Split
  !
  interface PAR_COMM_SPLIT
     module procedure PAR_COMM_SPLIT4,&
          &           PAR_COMM_SPLIT8,&
          &           PAR_COMM_SPLIT_IP
  end interface PAR_COMM_SPLIT
  !
  ! Free
  !
  interface PAR_COMM_FREE
     module procedure PAR_COMM_FREE4,&
          &           PAR_COMM_FREE8
  end interface PAR_COMM_FREE
    
  public :: PAR_DEFINE_COMMUNICATOR            ! Define the communicator according to some keywords
  public :: PAR_COMM_RANK_AND_SIZE             ! Give rank (and size) of a communicator
  public :: PAR_INIT                           ! Initialize MPI
  public :: PAR_LENGTH_INTEGER                 ! Length of integers
  public :: PAR_COMM_SPLIT                     ! Split a communicator
  public :: PAR_COMM_FREE                      ! Free MPI communicator
  public :: PAR_BARRIER                        ! Barrier
  public :: PAR_MPI_ERROR_TO_MESSAGE           ! Transform an MPI error code into a string
  public :: PAR_MPI_RUNEND                     ! End with an MPI message
  public :: PAR_COMM_SET_ERRHANDLER            ! Set error handler
  
  public :: PAR_IMASTER_IN_COMMUNICATOR        ! If I am a master of a given communicator
  public :: PAR_COMM_SPLIT_SIMVIZ              ! communicator splitting for viz ans im processes

contains
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-12-29
  !> @brief   Define communcator and communication arrays
  !> @details Define the communicator according to a keyword
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    implicit none

    character(*),        optional,          intent(in)  :: wherein
    integer(4),                             intent(out) :: PAR_COMM_TO_USE
    type(comm_data_par), optional, pointer, intent(inout) :: commu
    integer(ip)                                         :: icolo,jcolo

    !if( IPARALL ) then
    if( present(wherein) ) then
       if( trim(wherein) == 'IN THE UNIVERSE' ) then
          !
          ! In the universe
          !
#ifndef MPI_OFF
          PAR_COMM_TO_USE = MPI_COMM_WORLD
#endif
       else if( trim(wherein) == 'IN THE WORLD' ) then
          !
          ! In the world
          !
          PAR_COMM_TO_USE = int(PAR_COMM_WORLD,4_ip)   ! Alya world
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) )   commu => commd

       else if( trim(wherein) == 'IN MY CODE' ) then
          !
          ! In my code
          !
          !icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,4_ip)
          !PAR_COMM_COLOR(icolo,icolo)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_MY_CODE_ARRAY(1)

       else if( trim(wherein) == 'IN MY CODE WITHOUT MASTER' ) then
          !
          ! In my code without master
          !
          PAR_COMM_TO_USE = PAR_COMM_MY_CODE_WM4
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_MY_CODE_ARRAY(1)

       else if( trim(wherein) == 'IN MY ZONE' .or. trim(wherein) == 'IN CURRENT ZONE' ) then
          !
          ! In my current zone
          !
          icolo = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,icolo),4)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)

       else if( trim(wherein) == 'IN MY SUBD' ) then
          !
          ! In my current subd
          !
          icolo = par_code_zone_subd_to_color(current_code,0_ip,current_subd)
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,icolo),4)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)

       else if( trim(wherein) == 'IN CURRENT COUPLING' ) then
          !
          ! In my current coupling
          !
          icolo = color_target
          jcolo = color_source
          ! REVISAR IN CURRENT COUPLING
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,jcolo),4)
          if( present(commu) ) then
             call runend('PAR_DEFINE_COMMUNICATOR: WRONG OPTION 1')
          end if

       else if( trim(wherein) == 'IN CURRENT COLOR' .or. trim(wherein) == 'IN CURRENT TARGET COLOR' ) then
          !
          ! In my current target color
          !
          icolo = color_target
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,icolo),4)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)

       else if( trim(wherein) == 'IN CURRENT SOURCE COLOR' ) then
          !
          ! In my current source color
          !
          icolo = color_source
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,icolo),4)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)
          
       else if( trim(wherein) == 'IN CURRENT' ) then
          !
          ! Uses current communicator
          !
          PAR_COMM_TO_USE = int(PAR_COMM_CURRENT,4)

          if( present(commu) ) then
             call runend('PAR_DEFINE_COMMUNICATOR: WRONG OPTION 2')
          end if
       else if( trim(wherein) == 'IN SFC PARTITION' ) then

          PAR_COMM_TO_USE = int(PAR_COMM_SFC_WM,4)
       else if( trim(wherein) == 'IN SFC PARTITION WITH MASTER' ) then

          PAR_COMM_TO_USE = int(PAR_COMM_SFC,4)
       else if( trim(wherein) == 'IN MPIO' ) then

          PAR_COMM_TO_USE = int(PAR_COMM_MPIO_WM,4)
       else if( trim(wherein) == 'IN MPIO WITH MASTER' ) then

          PAR_COMM_TO_USE = int(PAR_COMM_MPIO,4)

       else

          call runend('PAR DEFINE COMMUNICATOR: INVALID COMMUNICATOR OPTION: '//trim(wherein))

       end if
    else
       PAR_COMM_TO_USE =  int(PAR_COMM_WORLD,4)
       commu            => commd
    end if
    !else
    !   PAR_COMM_TO_USE =  0
    !   commu            => commd
    !end if

  end subroutine PAR_DEFINE_COMMUNICATOR
  !-----------------------------------------------------------------------
  !
  !> @brief   Initialize MPI
  !> @details Initialize MPI
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_INIT()
    implicit none
    integer(4) :: istat4

#ifndef MPI_OFF
    istat4 = 0
    call MPI_Init(istat4)
    if( istat4 /= MPI_SUCCESS ) call runend('COULD NOT INITIALIZE MPI')
#endif
#ifdef EXTRAE
    call extrae_shutdown()
#endif

  end subroutine PAR_INIT

  !-----------------------------------------------------------------------
  !
  !> @brief   Define length of integers
  !> @details Define length of integers
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_LENGTH_INTEGER()
    implicit none

#ifndef MPI_OFF
    if( ip == 4 ) then
       PAR_INTEGER = MPI_INTEGER
    else
       PAR_INTEGER = MPI_INTEGER8
    end if
#endif

  end subroutine PAR_LENGTH_INTEGER

  !----------------------------------------------------------------------
  !
  ! SPLIT COMMUNICATOR
  ! IKEY should have the rank
  !
  !----------------------------------------------------------------------

  subroutine PAR_COMM_SPLIT4(icolor,PAR_COMM_FINAL,my_new_rank,wherein)
    implicit none 
    integer(4),   intent(in)           :: icolor
    integer(4),   intent(out)          :: PAR_COMM_FINAL
    integer(4),   intent(out)          :: my_new_rank
    character(*), intent(in), optional :: wherein
    integer(4)                         :: ikey4,istat4,jcolor4
    integer(4)                         :: PAR_COMM_INITIAL

#ifndef MPI_OFF
    istat4 = 0_4
    if( icolor == 0 ) then
       jcolor4 = MPI_UNDEFINED
    else
       jcolor4 = int(icolor,4)
    end if
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_INITIAL)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_INITIAL)
    end if
    call MPI_COMM_RANK(PAR_COMM_INITIAL,ikey4,istat4)

    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT4: ERROR 1')
    call MPI_COMM_SPLIT(PAR_COMM_INITIAL,jcolor4,ikey4,PAR_COMM_FINAL,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT4: ERROR 2')
    if( PAR_COMM_FINAL /= MPI_COMM_NULL ) then
       call MPI_COMM_RANK (PAR_COMM_FINAL,my_new_rank,istat4)
       if( istat4 /= MPI_SUCCESS ) call runend('PAR_COMM_SPLIT4: MPI ERROR 3')
    else
       my_new_rank = -1
    end if
#else
    PAR_COMM_FINAL = 0
#endif

  end subroutine PAR_COMM_SPLIT4

  subroutine PAR_COMM_SPLIT8(icolor,PAR_COMM_FINAL,my_new_rank,wherein)
    
    integer(8),   intent(in)           :: icolor
    integer(8),   intent(out)          :: PAR_COMM_FINAL
    integer(8),   intent(out)          :: my_new_rank
    character(*), intent(in), optional :: wherein
    integer(4)                         :: my_new_rank4
    integer(4)                         :: ikey4,istat4,jcolor4
    integer(4)                         :: PAR_COMM_INITIAL4
    integer(4)                         :: PAR_COMM_FINAL4

#ifndef MPI_OFF
    istat4 = 0_4
    if( icolor == 0 ) then
       jcolor4 = MPI_UNDEFINED
    else
       jcolor4 = int(icolor,4)
    end if

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_INITIAL4)
    call MPI_COMM_RANK (PAR_COMM_INITIAL4,ikey4,istat4)

    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT8: ERROR 1')
    call MPI_COMM_SPLIT(PAR_COMM_INITIAL4,jcolor4,ikey4,PAR_COMM_FINAL4,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT8: ERROR 2')
    if( PAR_COMM_FINAL4 /= MPI_COMM_NULL ) then
       call MPI_COMM_RANK (PAR_COMM_FINAL4,my_new_rank4,istat4)
       if( istat4 /= MPI_SUCCESS ) call runend('PAR_COMM_SPLIT8: MPI ERROR 3')
       PAR_COMM_FINAL = int(PAR_COMM_FINAL4,8)
       my_new_rank    = int(my_new_rank4,8)
    else
       my_new_rank    = -1
    end if
#else
    PAR_COMM_FINAL = 0
#endif

  end subroutine PAR_COMM_SPLIT8

  subroutine PAR_COMM_SPLIT_IP(icolor,PAR_COMM_FINAL,RANK_FINAL,PAR_COMM_INITIAL,RANK_INITIAL)
   
    integer(ip),   intent(in)           :: icolor
    integer(ip),   intent(out)          :: PAR_COMM_FINAL
    integer(ip),   intent(out)          :: RANK_FINAL
    integer(ip),   intent(in)           :: PAR_COMM_INITIAL
    integer(ip),   intent(in)           :: RANK_INITIAL
    integer(4)                          :: ikey4,istat4,jcolor4

#ifndef MPI_OFF
    istat4 = 0_4
    if( icolor == 0 ) then
       jcolor4 = MPI_UNDEFINED
    else
       jcolor4 = int(icolor,4)
    end if
    ikey4 = int(RANK_INITIAL,4)
    call MPI_COMM_SPLIT(PAR_COMM_INITIAL,jcolor4,ikey4,PAR_COMM_FINAL,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT_IP: ERROR 1')
    if( PAR_COMM_FINAL /= MPI_COMM_NULL ) then
       call MPI_COMM_RANK (PAR_COMM_FINAL,RANK_FINAL,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_SPLIT_IP: ERROR 2')
    else
       RANK_FINAL = -1
    end if
#else
    PAR_COMM_FINAL = 0
#endif

  end subroutine PAR_COMM_SPLIT_IP

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-06-11
  !> @brief   Free MPI communicator
  !> @details Free MPI communicator
  !> 
  !-----------------------------------------------------------------------
  
  subroutine PAR_COMM_FREE4(PAR_COMM_FINAL)
    integer(4),   intent(inout) :: PAR_COMM_FINAL
    integer(4)                  :: istat4
#ifndef MPI_OFF
    if( PAR_COMM_FINAL /= MPI_COMM_NULL ) call MPI_COMM_FREE(PAR_COMM_FINAL,istat4)
#endif    
  end subroutine PAR_COMM_FREE4
  
  subroutine PAR_COMM_FREE8(PAR_COMM_FINAL)
    integer(8),   intent(inout) :: PAR_COMM_FINAL
    integer(4)                  :: istat4
    integer(4)                  :: PAR_COMM_FINAL4
#ifndef MPI_OFF
    PAR_COMM_FINAL4 = int(PAR_COMM_FINAL,4)
    if( PAR_COMM_FINAL /= MPI_COMM_NULL ) call MPI_COMM_FREE(PAR_COMM_FINAL4,istat4)
#endif    
  end subroutine PAR_COMM_FREE8 

  !----------------------------------------------------------------------
  !
  ! RANK and SIZE of a communicator
  !
  !----------------------------------------------------------------------

  subroutine PAR_COMM_RANK_AND_SIZE_4(PAR_COMM_ORIGINAL,my_rank,comm_size)
    implicit none
    integer(4),            intent(in)  :: PAR_COMM_ORIGINAL  !< Communicator
    integer(4),            intent(out) :: my_rank            !< Rank in this communicator
    integer(4), optional,  intent(out) :: comm_size          !< Size of this communicator
    integer(4)                         :: istat4

#ifndef MPI_OFF
    istat4 = 0_4
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_4')
    end if
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_4')
#else
    my_rank = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_4

  subroutine PAR_COMM_RANK_AND_SIZE_8(PAR_COMM_ORIGINAL,my_rank,comm_size)
    implicit none
    integer(8),            intent(in)  :: PAR_COMM_ORIGINAL !< Communicator
    integer(8),            intent(out) :: my_rank           !< Rank in this communicator
    integer(8), optional,  intent(out) :: comm_size         !< Size of this communicator
    integer(4)                         :: istat4
    integer(4)                         :: my_rank4
    integer(4)                         :: comm_size4
    integer(4)                         :: PAR_COMM_ORIGINAL4

#ifndef MPI_OFF
    istat4 = 0_4
    PAR_COMM_ORIGINAL4 = int(PAR_COMM_ORIGINAL,4)
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_ORIGINAL4,comm_size4,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_8')
       comm_size = int(comm_size4,8)
    end if
    call MPI_Comm_rank(PAR_COMM_ORIGINAL4,my_rank4,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_8')
    my_rank = int(my_rank4,8)
#else
    my_rank = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_8

  subroutine PAR_COMM_RANK_AND_SIZE_4W(my_rank,comm_size,wherein)
    implicit none
    integer(4),             intent(out) :: my_rank           !< Rank in this communicator
    integer(4),             intent(out) :: comm_size         !< Size of this communicator
    character(*),           intent(in)  :: wherein             !< Wherein
    integer(4)                          :: istat4
    integer(4)                          :: PAR_COMM_ORIGINAL

#ifndef MPI_OFF
    istat4 = 0_4
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_ORIGINAL)
    call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_4W')
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_4W')
#else
    my_rank   = 0
    comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_4W

  subroutine PAR_COMM_RANK_AND_SIZE_41W(my_rank,wherein)
    implicit none
    integer(4),             intent(out) :: my_rank           !< Rank in this communicator
    character(*),           intent(in)  :: wherein             !< Wherein
    integer(4)                          :: istat4
    integer(4)                          :: PAR_COMM_ORIGINAL

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_ORIGINAL)
       call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_41W')
#endif
    else
       my_rank = -1
    end if

  end subroutine PAR_COMM_RANK_AND_SIZE_41W

  subroutine PAR_COMM_RANK_AND_SIZE_8W(my_rank,comm_size,wherein)
    implicit none
    integer(8),   intent(out) :: my_rank           !< Rank in this communicator
    integer(8),   intent(out) :: comm_size         !< Size of this communicator
    character(*), intent(in)  :: wherein           !< Wherein
    integer(4)                :: istat4
    integer(4)                :: PAR_COMM_ORIGINAL
    integer(4)                :: my_rank4
    integer(4)                :: comm_size4

#ifndef MPI_OFF
    istat4 = 0_4
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_ORIGINAL)
    call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size4,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_8W')
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank4,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_COMM_RANK_AND_SIZE_8W')
    my_rank = int(my_rank4,8)
    comm_size = int(comm_size4,8)
#else
    my_rank   = 0
    comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_8W

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   MPI Barrier
  !> @details MPI Barrier
  !
  !----------------------------------------------------------------------

  subroutine PAR_BARRIER(wherein,PAR_COMM_IN4)
    character(*), optional, intent(in) :: wherein
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(4)                         :: PAR_COMM_TO_USE4
    integer(4)                         :: istat4

#ifndef MPI_OFF
    if( IPARALL ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE4=PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE4)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE4)
          end if
       endif
       call MPI_Barrier( PAR_COMM_TO_USE4, istat4 )
    end if
#endif

  end subroutine PAR_BARRIER

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   Transform an MPI error code into a string
  !> @details Transform an MPI error code into a string
  !>
  !-----------------------------------------------------------------------

  function PAR_MPI_ERROR_TO_MESSAGE(istat4) result(message)

    integer(4),                    intent(in)  :: istat4   !< Error code
    character(len=:), allocatable              :: message   !< Message
    integer(4)                                 :: length4
    integer(4)                                 :: temp4

#ifndef MPI_OFF
    character(len=MPI_MAX_ERROR_STRING) :: message_mpi
    if( istat4 /= MPI_SUCCESS ) then
       call MPI_Error_string(istat4,message_mpi,length4,temp4)
       message = message_mpi(1:length4)
    end if
#endif

  end function PAR_MPI_ERROR_TO_MESSAGE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/04/2014
  !> @brief   Define if I am the master in this comunicatior
  !> @details Answer is TRUE if my rank is zero in communicator
  !
  !----------------------------------------------------------------------

  function PAR_IMASTER_IN_COMMUNICATOR(PAR_COMM_TO_USE)
    integer(4),  intent(in) :: PAR_COMM_TO_USE              !< Communicator
    integer(4)              :: my_rank
    logical(lg)             :: PAR_IMASTER_IN_COMMUNICATOR

    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank)
    if( my_rank == 0_4 ) then
       PAR_IMASTER_IN_COMMUNICATOR = .true.
    else
       PAR_IMASTER_IN_COMMUNICATOR = .false.
    end if

  end function PAR_IMASTER_IN_COMMUNICATOR

  subroutine PAR_COMM_SPLIT_SIMVIZ(PAR_COMM_SIM,PAR_COMM_VIZ,PAR_COMM_SIMVIZ,rank,mpisplit)

    implicit none
    integer(4) :: PAR_COMM_SIM, PAR_COMM_VIZ, PAR_COMM_SIMVIZ, rank, mpisplit
    integer(4) :: vizcolor,simcolor,istat4

#ifndef MPI_OFF

    if ( mod(rank+1,mpisplit) .ne. 0_ip ) then
        simcolor = 1
        vizcolor = MPI_UNDEFINED
     else
        simcolor = MPI_UNDEFINED
        vizcolor = 1
     end if

     call MPI_COMM_SPLIT(PAR_COMM_SIMVIZ,simcolor,rank,PAR_COMM_SIM,istat4)
     if( istat4 /= MPI_SUCCESS ) call runend('PAR_COMM_SPLIT_simviz: MPI ERROR 2')
     call MPI_COMM_SPLIT(PAR_COMM_SIMVIZ,vizcolor,rank,PAR_COMM_VIZ,istat4)
     if( istat4 /= MPI_SUCCESS ) call runend('PAR_COMM_SPLIT_simviz: MPI ERROR 2')

#endif

   end subroutine PAR_COMM_SPLIT_SIMVIZ
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   End
  !> @details End with an MPI error message
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_MPI_RUNEND(istat4,vacal) 

    integer(4),                 intent(in) :: istat4    !< Error code
    character(len=*), optional, intent(in) :: vacal     !< Caller
    
#ifndef MPI_OFF
    if( istat4 /= MPI_SUCCESS ) then
       if( present(vacal) ) then
          call runend(trim(vacal)//': '//PAR_MPI_ERROR_TO_MESSAGE(istat4))
       else
          call runend('MPI ERROR: '//PAR_MPI_ERROR_TO_MESSAGE(istat4))
       end if
    end if
#endif

  end subroutine PAR_MPI_RUNEND

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   End
  !> @details End with an MPI error message
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_COMM_SET_ERRHANDLER()

    integer(4) :: istat4
    
#ifndef MPI_OFF
      call MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN,istat4)
#endif

    end subroutine PAR_COMM_SET_ERRHANDLER
    
end module mod_communications_tools
!> @}
