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

module mod_communications_global

  use def_communications
  use mod_communications_tools
  implicit none  

  private
  
  integer(4)            :: ireq41(1)
  real(rp), allocatable :: xx_non_blocking(:)  
  !
  ! All reduce operations:
  !
  ! PAR_SUM
  ! PAR_MAX
  ! PAR_MIN
  !
  interface PAR_SUM
     module procedure PAR_SUM_IP_0,PAR_SUM_IP_s,PAR_SUM_IP_1,PAR_SUM_IP_2,PAR_SUM_IP_3,PAR_SUM_IP_03,&
          &           PAR_SUM_RP_0,PAR_SUM_RP_0b,PAR_SUM_RP_s,PAR_SUM_RP_1,PAR_SUM_RP_2,PAR_SUM_RP_3,&
          &           PAR_SUM_RP_02,&
          &           PAR_SUM_CX_0,PAR_SUM_CX_s,PAR_SUM_CX_1,PAR_SUM_CX_2
  end interface PAR_SUM
  interface PAR_MAX
     module procedure PAR_MAX_4_s,PAR_MAX_8_s,PAR_MAX_IP_0,PAR_MAX_IP_1,PAR_MAX_IP_2,&
          &           PAR_MAX_RP_s,PAR_MAX_RP_0,PAR_MAX_RP_1,PAR_MAX_RP_2,PAR_MAX_RP_3,&
          &           PAR_MAX_CX_s,PAR_MAX_CX_0,PAR_MAX_CX_1,PAR_MAX_CX_2
  end interface PAR_MAX
  interface PAR_MIN
     module procedure PAR_MIN_4_s ,PAR_MIN_8_s ,PAR_MIN_IP_0,PAR_MIN_IP_1,PAR_MIN_IP_2,&
          &           PAR_MIN_RP_s,PAR_MIN_RP_0,PAR_MIN_RP_1,PAR_MIN_RP_2,&
          &           PAR_MIN_CX_s,PAR_MIN_CX_0,PAR_MIN_CX_1,PAR_MIN_CX_2
  end interface PAR_MIN
  interface PAR_OR
     module procedure PAR_OR_LG_s
  end interface PAR_OR
  interface PAR_AVERAGE
     module procedure PAR_AVERAGE_IP_s,PAR_AVERAGE_IP_0,&
          &           PAR_AVERAGE_RP_s,PAR_AVERAGE_RP_0
  end interface PAR_AVERAGE
  interface PAR_LOAD_BALANCE
     module procedure PAR_LOAD_BALANCE_RP_0,PAR_LOAD_BALANCE_RP_s
  end interface PAR_LOAD_BALANCE
  !
  ! Reduction operations in the Alya world (PAR_COMM_WORLD) including masters
  !
  interface PAR_MAX_ALL
     module procedure PAR_MAX_ALL_IP_s, PAR_MAX_ALL_IP_1
  end interface PAR_MAX_ALL
  interface PAR_SUM_ALL
     module procedure PAR_SUM_ALL_IP_1, PAR_SUM_ALL_IP_3
  end interface PAR_SUM_ALL
  !
  ! All to all
  !
  interface PAR_ALLTOALL
     module procedure PAR_ALLTOALL_IP_0,&
          &           PAR_ALLTOALL_IP_14,&
          &           PAR_ALLTOALL_IP_18
  end interface PAR_ALLTOALL
  !
  ! ALL_TO_ALLV
  !
  interface PAR_ALLTOALLV
     module procedure PAR_ALLTOALLV_IP_1,&
          &           PAR_ALLTOALLV_IP_2,&
          &           PAR_ALLTOALLV_RP_1,&
          &           PAR_ALLTOALLV_RP_2
  end interface PAR_ALLTOALLV
  !
  ! BORADCAST
  !
  interface PAR_BROADCAST_IP
     module procedure PAR_BROADCAST_I4,PAR_BROADCAST_I8
  end interface PAR_BROADCAST_IP

  interface PAR_BROADCAST
     module procedure PAR_BROADCAST_IP_04,PAR_BROADCAST_IP_08,&
          &           PAR_BROADCAST_I4_s,PAR_BROADCAST_I8_s,&
          &           PAR_BROADCAST_I4_1,PAR_BROADCAST_I8_1,&
          &           PAR_BROADCAST_R8_s,PAR_BROADCAST_R8_0,PAR_BROADCAST_R8_1,PAR_BROADCAST_R8_2,&
          &           PAR_BROADCAST_R4_s,PAR_BROADCAST_R4_0,PAR_BROADCAST_R4_1,PAR_BROADCAST_R4_2,&
          &           PAR_BROADCAST_R4_02,PAR_BROADCAST_R8_02,&
          &           PAR_BROADCAST_R16_s,PAR_BROADCAST_R16_0,PAR_BROADCAST_R16_1,&
          &           PAR_BROADCAST_LG_s,PAR_BROADCAST_LG_1,&
          &           PAR_BROADCAST_CH,PAR_BROADCAST_CH_1
  end interface PAR_BROADCAST
  !
  ! Operation on arrays between all partitions
  !
  interface PAR_ALL_TO_ALL_ARRAY_OPERATION
     module procedure PAR_ALL_TO_ALL_ARRAY_OPERATION_IP
  end interface PAR_ALL_TO_ALL_ARRAY_OPERATION
  !
  ! GATHER
  !
  interface PAR_GATHER
     module procedure PAR_GATHER_CHARACTER,&
          &           PAR_GATHER_IP_s4,PAR_GATHER_IP_s8,PAR_GATHER_IP_s48,PAR_GATHER_IP_14,PAR_GATHER_IP_18,&
          &           PAR_GATHER_IP_12,PAR_GATHER_IP_23,&
          &           PAR_GATHER_RP_s,PAR_GATHER_RP_1,PAR_GATHER_RP_12
  end interface PAR_GATHER
  !
  ! GATHERV
  !
  interface PAR_GATHERV
     module procedure PAR_GATHERV_RP_1,PAR_GATHERV_RP_21,PAR_GATHERV_RP_22,PAR_GATHERV_RP_33,PAR_GATHERV_RP_0,&
          &           PAR_GATHERV_IP_1,PAR_GATHERV_IP_21, PAR_GATHERV_IP_22,&
                      PAR_GATHERV_RP_21_SEND,PAR_GATHERV_RP_22_SEND,&
          &           PAR_GATHERV_IP_1_SEND,PAR_GATHERV_IP_21_SEND,PAR_GATHERV_IP_22_SEND
  end interface PAR_GATHERV
  !
  ! ALLGATHER
  !
  interface PAR_ALLGATHER
     module procedure PAR_ALLGATHER_s4,PAR_ALLGATHER_s8,      &
          &           PAR_ALLGATHER_IP_14,PAR_ALLGATHER_IP_18,&
          &           PAR_ALLGATHER_RP_0,PAR_ALLGATHER_RP_2,  &
          &           PAR_ALLGATHER_RP,                       &
          &           PAR_ALLGATHER_RP_02,PAR_ALLGATHER_LG
  end interface PAR_ALLGATHER  
  !
  ! ALLGATHERV
  !
  interface PAR_ALLGATHERV
     module procedure PAR_ALLGATHERV_IP4  ,PAR_ALLGATHERV_IP8  ,&
          &           PAR_ALLGATHERV_IP4_2,&
          &           PAR_ALLGATHERV_RP_1,&
          &           PAR_ALLGATHERV_RP_18,&
          &           PAR_ALLGATHERV_RP_2,&
          &           PAR_ALLGATHERV_RP_3
  end interface PAR_ALLGATHERV
  !
  ! SCATTER
  !
  interface PAR_SCATTER
     module procedure PAR_SCATTER_IP_s
  end interface PAR_SCATTER
  !
  ! SCATTERV
  !
  interface PAR_SCATTERV
     module procedure PAR_SCATTERV_IP_1,PAR_SCATTERV_IP_2,&
          &           PAR_SCATTERV_RP_1,PAR_SCATTERV_RP_2,PAR_SCATTERV_RP_0,&
          &           PAR_SCATTERV_IP_1_RCV,PAR_SCATTERV_IP_2_RCV,&
          &           PAR_SCATTERV_R4_2_RCV, PAR_SCATTERV_R8_2_RCV, PAR_SCATTERV_R16_2_RCV
  end interface PAR_SCATTERV
  
  public :: PAR_SUM                            ! AllReduce SUM
  public :: PAR_MAX                            ! AllReduce MAX
  public :: PAR_MIN                            ! AllReduce MIN
  public :: PAR_OR                             ! AllReduce OR
  public :: PAR_AVERAGE                        ! Average value
  public :: PAR_LOAD_BALANCE                   ! Load balance: average / max
  public :: PAR_SUM_ALL                        ! AllReduce SUM in the Alya world (PAR_COMM_WORLD)
  public :: PAR_MAX_ALL                        ! AllReduce MAX in the Alya world (PAR_COMM_WORLD)
  public :: PAR_ALLTOALL                       ! Alltoall
  public :: PAR_ALLTOALLV                      ! All to all V
  public :: PAR_ALLTOALLV_RP_1
  public :: PAR_ALLTOALLV_RP_2
  public :: PAR_ALLTOALLV_RP_4
  public :: PAR_ALLTOALLV_RP_3
  public :: PAR_ALLTOALLV_IP_1
  public :: PAR_ALLTOALLV_IP_2
  public :: PAR_BROADCAST                      ! Broadcast
  public :: PAR_ALL_TO_ALL_ARRAY_OPERATION     ! Array operations between all partitions of communicators
  public :: PAR_GATHER                         ! Gather of integers, reals, characters
  public :: PAR_GATHERV                        ! Gather of integers, reals, characters
  public :: PAR_ALLGATHERV                     ! All Gatherv
  public :: PAR_ALLGATHER                      ! All Gather
  public :: PAR_SCATTER                        ! Scatter of integers, reals, characters
  public :: PAR_SCATTERV                       ! Scatter of integers, reals, characters
  public :: PAR_WAITALL                        ! Waitall for non-blocking reductions
  
contains

  !----------------------------------------------------------------------
  !
  ! SUM FOR REALS
  !
  !----------------------------------------------------------------------

  subroutine PAR_SUM_RP_s(xx,wherein,commu)
    implicit none
    real(rp),                      intent(inout) :: xx
    character(*),        optional                :: wherein
    type(comm_data_par), optional, intent(in)    :: commu
    integer(ip)                                  :: nsize
    integer(4)                                   :: istat4,nsize4
    integer(4)                                   :: PAR_COMM_TO_USE
    real(rp)                                     :: yy
    character(30)                                :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(commu) ) then
          PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( PAR_COMM_TO_USE /= -1 ) then
          nsize  = 1
          nsize4 = int(nsize,4)
          if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
             yy = 0.0_rp
          else
             yy = xx
          end if
          call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,istat4)
          if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_s') 
       end if
    end if
#endif

  end subroutine PAR_SUM_RP_s

  subroutine PAR_SUM_RP_0(n,xx,wherein,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)    :: n
    real(rp),     intent(inout) :: xx(n)
    character(*), optional      :: wherein
    character(*), optional      :: wsynch
    integer(4),   optional      :: PAR_COMM_IN4
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    logical(lg)                 :: asynch
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE=PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       endif
       asynch = .false.
       if( present(wsynch) ) then
          if( trim(wsynch) == 'NON BLOCKING' ) asynch = .true.
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(xx_non_blocking(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             xx_non_blocking(ii) = 0.0_rp
          end do
       else
          do ii = 1,nsize
             xx_non_blocking(ii) = xx(ii)
          end do
       end if
       if( asynch ) then
#ifdef MPI3
          call MPI_IAllReduce(xx_non_blocking,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,ireq41(1),istat4)
#else
          call runend('THUS FUNCTIONALITY REQUIRES MPI3: COMPILE ALYA WITH -DMPI3')
#endif
       else
          call MPI_AllReduce(xx_non_blocking,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,istat4)
          deallocate(xx_non_blocking)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_0') 
    end if
#endif

  end subroutine PAR_SUM_RP_0

  subroutine PAR_SUM_RP_0b(xx,PAR_COMM_IN4)
    implicit none
    real(rp),     intent(inout) :: xx
    integer(4),   intent(in)    :: PAR_COMM_IN4
    real(rp)                    :: aa(1)

    aa(1)=xx
    call PAR_SUM_RP_0(1_ip,aa,PAR_COMM_IN4=PAR_COMM_IN4)
    xx=aa(1)

  end subroutine PAR_SUM_RP_0b

  subroutine PAR_SUM_RP_02(n1,n2,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n1
    integer(ip),  intent(in)    :: n2
    real(rp),     intent(inout) :: xx(n1,n2)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,jj,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    real(rp),     allocatable   :: yy(:,:)
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    nsize  = n1*n2
    if( IPARALL .and. nsize > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize4 = int(nsize,4)
       allocate(yy(n1,n2))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = 1,n2
             do ii = 1,n1
                yy(ii,jj) = 0.0_rp
             end do
          end do
       else
          do jj = 1,n2
             do ii = 1,n1
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_02') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_RP_02

  subroutine PAR_SUM_RP_1(xx,wherein,wsynch,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:)
    character(*), optional,intent(in)    :: wherein
    character(*), optional,intent(in)    :: wsynch
    integer(4),   optional, intent(in)   :: PAR_COMM_IN4
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    logical(lg)                          :: asynch
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       endif
       asynch = .false.
       if( present(wsynch) ) then
          if( trim(wsynch) == 'NON BLOCKING' ) asynch = .true.
       end if
       nsize  = memory_size(xx)
       nsize4 = int(nsize,4)
       allocate( xx_non_blocking(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(xx_non_blocking,1),ubound(xx_non_blocking,1)
             xx_non_blocking(ii) = 0.0_rp
          end do
       else
          do ii = lbound(xx_non_blocking,1),ubound(xx_non_blocking,1)
             xx_non_blocking(ii) = xx(ii)
          end do
       end if
      if( asynch ) then
#ifdef MPI3
          call MPI_IAllReduce(xx_non_blocking,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,ireq41(1),istat4)
#else
          call runend('THUS FUNCTIONALITY REQUIRES MPI3: COMPILE ALYA WITH -DMPI3')
#endif
       else
          call MPI_AllReduce(xx_non_blocking,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,istat4)
          deallocate( xx_non_blocking )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_1') 
    end if
#endif

  end subroutine PAR_SUM_RP_1

  subroutine PAR_SUM_RP_2(xx,wherein)
    implicit none
    real(rp),     pointer, intent(inout)         :: xx(:,:)
    character(*),          optional              :: wherein
    integer(ip)                                  :: ii,jj,nsize1,nsize2,nsize
    integer(4)                                   :: istat4,nsize4
    integer(4)                                   :: PAR_COMM_TO_USE
    real(rp),     pointer                        :: yy(:,:)
    character(30)                                :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = memory_size(xx,1_ip)
       nsize2 = memory_size(xx,2_ip)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = 0.0_rp
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_2') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_RP_2

  subroutine PAR_SUM_RP_3(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer, intent(inout)         :: xx(:,:,:)
    character(*),          optional              :: wherein
    character(*), optional                       :: who
    integer(4),   optional, intent(in)           :: PAR_COMM_IN4
    integer(ip)                                  :: ii,jj
    integer(4)                                   :: istat4,nsize4
    integer(4)                                   :: PAR_COMM_TO_USE
    real(rp),     pointer                        :: yy(:,:,:)
    logical(lg)                                  :: not_include_master
    character(30)                                :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize4 = int(size(xx),4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2),lbound(xx,3):ubound(xx,3)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0.0_rp
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_3') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_RP_3

  !----------------------------------------------------------------------
  !
  ! SUM FOR COMPLEX
  !
  !----------------------------------------------------------------------

  subroutine PAR_SUM_CX_s(xx,wherein)
    implicit none
    complex(rp),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp)                 :: yy
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          yy = 0.0_rp
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_CX_s') 
    end if
#endif

  end subroutine PAR_SUM_CX_s

  subroutine PAR_SUM_CX_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    complex(rp),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp),  allocatable   :: yy(:)
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) = 0.0_rp
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_CX_0') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_CX_0

  subroutine PAR_SUM_CX_1(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( associated(xx) ) then
          nsize  = size(xx,1)
       else
          nsize = 0
       end if
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = 0.0_rp
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_CX_1') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_CX_1

  subroutine PAR_SUM_CX_2(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:,:)
    character(*), optional               :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:,:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( associated(xx) ) then
          nsize1 = size(xx,1)
          nsize2 = size(xx,2)
       else
          nsize1 = 0
          nsize2 = 0         
       end if
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = (0.0_rp,0.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_CX_2')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_CX_2

  !----------------------------------------------------------------------
  !
  ! AVERAGE FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_AVERAGE_IP_s(xx,wherein)
    implicit none
    integer(ip),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(xx,wherein)
    else
       call PAR_SUM(xx)
    end if
    if( IPARALL ) xx = xx / (max(1_ip,comm_size4-1_ip))

  end subroutine PAR_AVERAGE_IP_s

  subroutine PAR_AVERAGE_IP_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    integer(ip),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(n,xx,wherein)
    else
       call PAR_SUM(n,xx)
    end if
    if( IPARALL ) xx(1:n) = xx(1:n) / (max(1_ip,comm_size4-1_ip))

  end subroutine PAR_AVERAGE_IP_0

  !----------------------------------------------------------------------
  !
  ! AVERAGE FOR REAL
  !
  !----------------------------------------------------------------------

  subroutine PAR_AVERAGE_RP_s(xx,wherein)
    implicit none
    real(rp),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(xx,wherein)
    else
       call PAR_SUM(xx)
    end if
    if( IPARALL ) xx = xx / real(max(1_ip,comm_size4-1_ip),rp)

  end subroutine PAR_AVERAGE_RP_s

  subroutine PAR_AVERAGE_RP_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    real(rp),     intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(n,xx,wherein)
    else
       call PAR_SUM(n,xx)
    end if
    if( IPARALL ) xx(1:n) = xx(1:n) / real(max(1_ip,comm_size4-1_ip),rp)

  end subroutine PAR_AVERAGE_RP_0

 !----------------------------------------------------------------------
  !
  ! AVERAGE FOR REAL
  !
  !----------------------------------------------------------------------

  subroutine PAR_LOAD_BALANCE_RP_s(xx,wherein)
    implicit none
    real(rp),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4
    real(rp)                    :: xx_max

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    call PAR_SUM(xx,wherein)
    xx_max = xx
    if( present(wherein) ) then
       call PAR_MAX(xx_max,wherein)
    else
       call PAR_MAX(xx_max)
    end if
    if( IPARALL ) xx = xx / real(max(1_ip,comm_size4-1_ip),rp)
    xx = xx / (xx_max + epsilon(1.0_rp))

  end subroutine PAR_LOAD_BALANCE_RP_s

  subroutine PAR_LOAD_BALANCE_RP_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    real(rp),     intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4
    real(rp)                    :: xx_max(n)

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(n,xx,wherein)
    else
       call PAR_SUM(n,xx)
    end if
    xx_max = xx
    if( IPARALL ) xx(1:n) = xx(1:n) / real(max(1_ip,comm_size4-1_ip),rp)
     xx(1:n) = xx(1:n) / (xx_max(1:n) + epsilon(1.0_rp))

  end subroutine PAR_LOAD_BALANCE_RP_0

  !----------------------------------------------------------------------
  !
  ! SUM FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_SUM_IP_s(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),            intent(inout) :: xx
    character(*), optional                :: wherein
    character(*), optional                :: who
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(ip)                           :: nsize
    integer(4)                            :: istat4,nsize4
    integer(4)                            :: PAR_COMM_TO_USE
    integer(ip)                           :: yy
    logical(lg)                           :: not_include_master
    character(30)                         :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0_ip
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_IP_s') 
    end if
#endif

  end subroutine PAR_SUM_IP_s

  subroutine PAR_SUM_IP_0(n,xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)            :: n
    integer(ip),  intent(inout)         :: xx(n)
    character(*), intent(in),  optional :: wherein
    character(*), intent(in),  optional :: who
    integer(4),   intent(in),  optional :: PAR_COMM_IN4
    integer(ip)                         :: ii,nsize
    integer(4)                          :: istat4,nsize4
    integer(4)                          :: PAR_COMM_TO_USE
    integer(ip),  allocatable           :: yy(:)
    logical(lg)                         :: not_include_master
    character(30)                       :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0_ip
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_IP_0') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_IP_0

  subroutine PAR_SUM_IP_03(n1,n2,n3,xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)           :: n1,n2,n3
    integer(ip),  intent(inout)        :: xx(n1,n2,n3)
    character(*), optional             :: wherein
    character(*), optional             :: who
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(ip)                        :: ii
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    integer(ip),  allocatable          :: yy(:,:,:)
    logical(lg)                        :: not_include_master
    character(30)                      :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n1*n2*n3 > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize4 = int(n1*n2*n3,4)
       allocate(yy(n1,n2,n3))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0_ip
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_IP_03') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_IP_03

  subroutine PAR_SUM_IP_1(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),     pointer,  intent(inout) :: xx(:)
    character(*),    optional                :: wherein
    character(*),    optional                :: who
    integer(4),      optional, intent(in)    :: PAR_COMM_IN4
    integer(ip)                              :: ii,nsize
    integer(4)                               :: istat4,nsize4
    integer(4)                               :: PAR_COMM_TO_USE
    integer(ip),     pointer                 :: yy(:)
    logical(lg)                              :: not_include_master
    character(30)                            :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize  = memory_size(xx,1_ip)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = 0_ip
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_IP_1') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_IP_1

  subroutine PAR_SUM_IP_2(xx,wherein,who)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:,:)
    character(*),    optional               :: wherein
    character(*),    optional               :: who
    integer(ip)                             :: ii,jj,nsize1,nsize2,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:,:)
    logical(lg)                             :: not_include_master
    character(30)                           :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if

       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize1 = memory_size(xx,1_ip)
       nsize2 = memory_size(xx,2_ip)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = 0_ip
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_IP_2') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_IP_2

  subroutine PAR_SUM_IP_3(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),    pointer, intent(inout) :: xx(:,:,:)
    character(*),   optional               :: wherein
    character(*),   optional               :: who
    integer(4),     optional, intent(in)   :: PAR_COMM_IN4
    integer(4)                             :: istat4,nsize4
    integer(4)                             :: PAR_COMM_TO_USE
    integer(ip),    pointer                :: yy(:,:,:)
    logical(lg)                            :: not_include_master
    character(30)                          :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize4 = int(size(xx),4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2),lbound(xx,3):ubound(xx,3)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_IP_3') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_IP_3

  !----------------------------------------------------------------------
  !
  ! MAX FOR REALS
  !
  !----------------------------------------------------------------------

  subroutine PAR_MAX_RP_s(xx,wherein,who,rank_max_owner,PAR_COMM_IN4)
    implicit none
    real(rp),     intent(inout) :: xx
    character(*), optional      :: wherein
    character(*), optional      :: who
    integer(ip),  optional      :: rank_max_owner
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: rank_min
    integer(4)                  :: rank_max_owner4
    real(rp)                    :: yy
    logical(lg)                 :: not_include_master
    character(30)               :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
        if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             my_wherein = 'IN MY CODE'
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = -huge(1.0_rp)
       else
          yy = xx
       end if

       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_s') 
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( abs(yy-xx) < epsilon(1.0_rp) ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,ip)
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MAX_RP_s

  subroutine PAR_MAX_RP_0(n,xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)           :: n
    real(rp),     intent(inout)        :: xx(n)
    character(*), optional, intent(in) :: wherein
    character(*), optional, intent(in) :: who
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(ip)                        :: ii,nsize
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    real(rp),     allocatable          :: yy(:)
    logical(lg)                        :: not_include_master
    character(30)                      :: my_wherein

    if( IPARALL .and. n > 0 ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
        if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do ii = 1,nsize
             yy(ii) = -huge(1.0_rp)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_0') 
       deallocate(yy)
#endif
    end if

  end subroutine PAR_MAX_RP_0

  subroutine PAR_MAX_RP_1(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:)
    character(*), optional, intent(in)   :: wherein
    character(*), optional, intent(in)   :: who
    integer(4),   optional, intent(in)   :: PAR_COMM_IN4
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    real(rp),     pointer                :: yy(:)
    logical(lg)                          :: not_include_master
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = memory_size(xx)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = -huge(1.0_rp)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_1') 
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_RP_1

  subroutine PAR_MAX_RP_2(xx,wherein)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:,:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    real(rp),     pointer                :: yy(:,:)
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = memory_size(xx,1_ip)
       nsize2 = memory_size(xx,2_ip)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = -huge(1.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_2') 
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_RP_2

  subroutine PAR_MAX_RP_3(xx,wherein)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:,:,:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,jj,kk,nsize1,nsize2,nsize3,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    real(rp),     pointer                :: yy(:,:,:)
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = memory_size(xx,1_ip)
       nsize2 = memory_size(xx,2_ip)
       nsize3 = memory_size(xx,3_ip)
       nsize  = nsize1*nsize2*nsize3
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2),lbound(xx,3):ubound(xx,3)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do kk = lbound(yy,3),ubound(yy,3)
             do jj = lbound(yy,2),ubound(yy,2)
                do ii = lbound(yy,1),ubound(yy,1)
                   yy(ii,jj,kk) = -huge(1.0_rp)
                end do
             end do
          end do
       else
          do kk = lbound(yy,3),ubound(yy,3)
             do jj = lbound(yy,2),ubound(yy,2)
                do ii = lbound(yy,1),ubound(yy,1)
                   yy(ii,jj,kk) = xx(ii,jj,kk)
                end do
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SUM_RP_3') 
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_RP_3

  !----------------------------------------------------------------------
  !
  ! MAX FOR COMPLEX
  !
  !----------------------------------------------------------------------

  subroutine PAR_MAX_CX_s(xx,wherein)
    implicit none
    complex(rp),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp)                 :: yy
    character(30)               :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          yy = -huge(1.0_rp)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_CX_s') 
#endif
    end if

  end subroutine PAR_MAX_CX_s

  subroutine PAR_MAX_CX_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    complex(rp),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp),  allocatable   :: yy(:)
    character(30)               :: my_wherein

    if( IPARALL .and. n > 0 ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) = -huge(1.0_rp)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_CX_0') 
       deallocate(yy)
#endif
    end if

  end subroutine PAR_MAX_CX_0

  subroutine PAR_MAX_CX_1(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:)
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( associated(xx) ) then
          nsize  = size(xx,1)
       else
          nsize = 0
       end if
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = -huge(1.0_rp)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_CX_1') 
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_CX_1

  subroutine PAR_MAX_CX_2(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:,:)
    character(*), optional               :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:,:)
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( associated(xx) ) then
          nsize1 = size(xx,1)
          nsize2 = size(xx,2)
       else
          nsize1 = 0
          nsize2 = 0
       end if
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = -huge(1.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_CX_2') 
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_CX_2

  !----------------------------------------------------------------------
  !
  ! MAX FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_MAX_4_s(xx,wherein,rank_max_owner,PAR_COMM_IN4,INCLUDE_ROOT)
    implicit none
    integer(4),   intent(inout)           :: xx
    character(*), intent(in),    optional :: wherein
    integer(4),   intent(inout), optional :: rank_max_owner
    integer(4),   intent(in),    optional :: PAR_COMM_IN4
    logical(lg),  intent(in),    optional :: INCLUDE_ROOT
    integer(ip)                           :: nsize
    integer(4)                            :: istat4,nsize4
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: rank_min
    integer(4)                            :: rank_max_owner4
    integer(4)                            :: yy
    logical(lg)                           :: if_include_root
    character(30)                         :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( present(INCLUDE_ROOT) ) then
          if_include_root = INCLUDE_ROOT
       else
          if_include_root = .false.
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. .not. if_include_root ) then
          yy = -huge(1_4)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_4_s') 

       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( xx == yy ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,ip)
       end if

#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MAX_4_s

  subroutine PAR_MAX_8_s(xx,wherein,rank_max_owner,PAR_COMM_IN4,INCLUDE_ROOT)
    implicit none
    integer(8),   intent(inout)           :: xx
    character(*), intent(in),    optional :: wherein
    integer(8),   intent(inout), optional :: rank_max_owner
    integer(4),   intent(in),    optional :: PAR_COMM_IN4
    logical(lg),  intent(in),    optional :: INCLUDE_ROOT
    integer(ip)                           :: nsize
    integer(4)                            :: istat4,nsize4
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: rank_min
    integer(4)                            :: rank_max_owner4
    integer(8)                            :: yy
    logical(lg)                           :: if_include_root
    character(30)                         :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( present(INCLUDE_ROOT) ) then
          if_include_root = INCLUDE_ROOT
       else
          if_include_root = .false.
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. .not. if_include_root ) then
          yy = -huge(1_8)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_INTEGER8,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_8_s') 
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( yy == xx ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,8)
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MAX_8_s

  subroutine PAR_MAX_IP_0(n,xx,wherein,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)    :: n
    integer(ip),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(4),   optional      :: PAR_COMM_IN4
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(ip),  allocatable   :: yy(:)
    character(30)               :: my_wherein

    if( IPARALL .and. n > 0 ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) = -huge(1_ip)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_IP_0') 
       deallocate(yy)
#endif
    end if

  end subroutine PAR_MAX_IP_0

  subroutine PAR_MAX_IP_1(xx,wherein)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:)
    character(*),    optional               :: wherein
    integer(ip)                             :: ii,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:)
    character(30)                           :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = memory_size(xx,1_ip)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = -huge(1_ip)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_IP_1') 
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_IP_1

  subroutine PAR_MAX_IP_2(xx,wherein)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:,:)
    character(*),    optional               :: wherein
    integer(ip)                             :: ii,jj,nsize1,nsize2,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:,:)
    character(30)                           :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = memory_size(xx,1_ip)
       nsize2 = memory_size(xx,2_ip)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = -huge(1_ip)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_IP_2') 
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_IP_2

  !----------------------------------------------------------------------
  !
  ! OR FOR REALS
  !
  !----------------------------------------------------------------------

  subroutine PAR_OR_LG_s(xx_lg,wherein)
    implicit none
    logical(lg),  intent(inout) :: xx_lg
    character(*), optional      :: wherein
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(ip)                 :: yy,xx
    character(30)               :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize4 = 1_4
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          yy =  0
       else
          if( xx_lg ) then
             yy = 1
          else
             yy = 0
          end if
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_OR_LS_s') 
       if( xx == 1 ) then
          xx_lg = .true.
       else
          xx_lg = .false.
       end if
#endif
    end if

  end subroutine PAR_OR_LG_s

  !----------------------------------------------------------------------
  !
  ! MIN FOR REALS
  !
  !----------------------------------------------------------------------

  subroutine PAR_MIN_RP_s(xx,wherein,who,rank_max_owner)
    implicit none
    real(rp),     intent(inout) :: xx
    character(*), optional      :: wherein
    character(*), optional      :: who
    integer(ip),  optional      :: rank_max_owner
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: rank_min
    integer(4)                  :: rank_max_owner4
    real(rp)                    :: yy
    logical(lg)                 :: not_include_master
    character(30)               :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy =  huge(1.0_rp)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_RP_s') 
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( abs(yy-xx) < epsilon(1.0_rp) ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,ip)
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MIN_RP_s

  subroutine PAR_MIN_RP_0(n,xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)           :: n
    real(rp),     intent(inout)        :: xx(n)
    character(*), optional, intent(in) :: wherein
    character(*), optional, intent(in) :: who
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(ip)                        :: ii,nsize
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    real(rp),     allocatable          :: yy(:)
    logical(lg)                        :: not_include_master
    character(30)                      :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL .and. n > 0 ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do ii = 1,nsize
             yy(ii) =  huge(1.0_rp)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_RP_0') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_MIN_RP_0

  subroutine PAR_MIN_RP_1(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer,  intent(inout) :: xx(:)
    character(*), optional, intent(in)    :: wherein
    character(*), optional, intent(in)    :: who
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(ip)                           :: ii,nsize
    integer(4)                            :: istat4,nsize4
    integer(4)                            :: PAR_COMM_TO_USE
    real(rp),     pointer                 :: yy(:)
    logical(lg)                           :: not_include_master
    character(30)                         :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize  = memory_size(xx)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) =  huge(1.0_rp)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_RP_1') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_MIN_RP_1

  subroutine PAR_MIN_RP_2(xx,wherein)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:,:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    real(rp),     pointer                :: yy(:,:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL ) then
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = memory_size(xx,1_ip)
       nsize2 = memory_size(xx,2_ip)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) =  huge(1.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_RP_2') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_MIN_RP_2

  !----------------------------------------------------------------------
  !
  ! MIN FOR COMPLEX
  !
  !----------------------------------------------------------------------

  subroutine PAR_MIN_CX_s(xx,wherein)
    implicit none
    complex(rp),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp)                 :: yy
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL ) then
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          yy = huge(1.0_rp)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_CX_s') 
    end if
#endif

  end subroutine PAR_MIN_CX_s

  subroutine PAR_MIN_CX_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    complex(rp),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp),  allocatable   :: yy(:)
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL .and. n > 0 ) then
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) = huge(1.0_rp)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_CX_0') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_MIN_CX_0

  subroutine PAR_MIN_CX_1(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( associated(xx) ) then
          nsize  = size(xx,1)
       else
          nsize = 0
       end if
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = huge(1.0_rp)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_CX_1') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_MIN_CX_1

  subroutine PAR_MIN_CX_2(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:,:)
    character(*), optional               :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:,:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( associated(xx) ) then
          nsize1 = size(xx,1)
          nsize2 = size(xx,2)
       else
          nsize1 = 0
          nsize2 = 0
       end if
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = huge(1.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_CX_2') 
       deallocate( yy )
    end if
#endif

  end subroutine PAR_MIN_CX_2

  !----------------------------------------------------------------------
  !
  ! MIN FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_MIN_4_s(xx,wherein,rank_max_owner,PAR_COMM_IN4,INCLUDE_ROOT)
    implicit none
    integer(4),   intent(inout)        :: xx
    character(*),             optional :: wherein
    integer(4),               optional :: rank_max_owner
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    logical(lg),  intent(in), optional :: INCLUDE_ROOT
    integer(ip)                        :: nsize
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    integer(4)                         :: yy
    integer(4)                         :: rank_min
    integer(4)                         :: rank_max_owner4
    logical(lg)                        :: if_include_root
    character(30)                      :: my_wherein

    istat4 = 0_4
    if( IPARALL ) then
#ifndef MPI_OFF
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( present(INCLUDE_ROOT) ) then
          if_include_root = INCLUDE_ROOT
       else
          if_include_root = .false.
       end if

       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. .not. if_include_root ) then
          yy =  huge(4_4)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_INTEGER,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_4_s') 
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( abs(yy-xx) == 0_4 ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = rank_max_owner4
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MIN_4_s

  subroutine PAR_MIN_8_s(xx,wherein,rank_max_owner,PAR_COMM_IN4,INCLUDE_ROOT)
    implicit none
    integer(8),   intent(inout)        :: xx
    character(*), optional             :: wherein
    integer(8),   optional             :: rank_max_owner
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    logical(lg),  intent(in), optional :: INCLUDE_ROOT
    integer(ip)                        :: nsize
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    integer(8)                         :: yy
    integer(4)                         :: rank_min
    integer(4)                         :: rank_max_owner4
    logical(lg)                        :: if_include_root
    character(30)                      :: my_wherein

    istat4 = 0_4
    if( IPARALL ) then
#ifndef MPI_OFF
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( present(INCLUDE_ROOT) ) then
          if_include_root = INCLUDE_ROOT
       else
          if_include_root = .false.
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. .not. if_include_root ) then
          yy =  huge(1_8)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_INTEGER8,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_8_s') 
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( abs(yy-xx) == 0_8 ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,8)
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MIN_8_s


  subroutine PAR_MIN_IP_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    integer(ip),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(ip),  allocatable   :: yy(:)
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) =  huge(1_ip)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_IP_0') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_MIN_IP_0

  subroutine PAR_MIN_IP_1(xx,wherein)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:)
    character(*),    optional               :: wherein
    integer(ip)                             :: ii,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:)
    character(30)                           :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = memory_size(xx,1_ip)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) =  huge(1_ip)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       deallocate( yy )
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_IP_1') 
    end if
#endif

  end subroutine PAR_MIN_IP_1

  subroutine PAR_MIN_IP_2(xx,wherein)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:,:)
    character(*),    optional               :: wherein
    integer(ip)                             :: ii,jj,nsize1,nsize2,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:,:)
    character(30)                           :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = memory_size(xx,1_ip)
       nsize2 = memory_size(xx,2_ip)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) =  huge(1_ip)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       deallocate( yy )
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_IP_2') 
    end if
#endif

  end subroutine PAR_MIN_IP_2
  !
  ! SUM for integers in the Alya world (PAR_COMM_WORLD) including masters
  !
  subroutine PAR_SUM_ALL_IP_1(xx)

    implicit none
    integer(ip),  pointer, intent(inout) :: xx(:)
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    integer(ip),  pointer                :: yy(:)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_TO_USE)
       nsize = memory_size(xx,1_ip)
       nsize4 = int(nsize,4)
       allocate(yy(lbound(xx,1):ubound(xx,1)))
       do ii =  lbound(yy,1), ubound(yy,1)
          yy(ii) = xx(ii)
       end do
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_IP_1') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_ALL_IP_1

  subroutine PAR_SUM_ALL_IP_3(xx)

    implicit none
    integer(ip),  pointer, intent(inout) :: xx(:,:,:)
    integer(ip)                          :: ii,jj,kk,nsize1,nsize2,nsize3,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    integer(ip),  pointer                :: yy(:,:,:)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_TO_USE)
       nsize1  = memory_size(xx,1_ip)
       nsize2  = memory_size(xx,2_ip)
       nsize3  = memory_size(xx,3_ip)
       nsize   = nsize1 * nsize2 * nsize3
       nsize4 = int(nsize,4)
       allocate(yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2),lbound(xx,3):ubound(xx,3)))
       do kk = lbound(yy,3), ubound(yy,3)
          do jj = lbound(yy,2), ubound(yy,2)
             do ii =  lbound(yy,1), ubound(yy,1)
                yy(ii,jj,kk) = xx(ii,jj,kk)
             end do
          end do
       end do
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MIN_IP_3') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_ALL_IP_3
  !
  ! MAX for integers in the Alya world (PAR_COMM_WORLD) including masters
  !
  subroutine PAR_MAX_ALL_IP_s(xx)

    implicit none
    integer(ip),  intent(inout) :: xx
    integer(4)                  :: istat4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(ip)                 :: yy

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_TO_USE)
       yy = xx
       call MPI_AllReduce(yy,xx,1_ip,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_ALL_IP_s') 
    end if
#endif
  end subroutine PAR_MAX_ALL_IP_s

  subroutine PAR_MAX_ALL_IP_1(xx)

    implicit none
    integer(ip),  pointer, intent(inout) :: xx(:)
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    integer(ip),  pointer                :: yy(:)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_TO_USE)
       nsize = memory_size(xx,1_ip)
       nsize4 = int(nsize,4)
       allocate(yy(lbound(xx,1):ubound(xx,1)))
       do ii =  lbound(yy,1), ubound(yy,1)
          yy(ii) = xx(ii)
       end do
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_MAX_ALL_IP_1') 
       deallocate(yy)
    end if
#endif

  end subroutine PAR_MAX_ALL_IP_1

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-06-04
  !> @brief   All to all
  !> @details MPI_ALLTOALL
  !>
  !-----------------------------------------------------------------------

!!$  subroutine PAR_ALLATOALL_IP_1(xx_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)
!!$
!!$    integer(ip),  pointer,  intent(in)  :: xx_send(:)
!!$    integer(ip),  pointer,  intent(inout) :: xx_recv(:)
!!$    character(*), optional, intent(in)  :: wherein
!!$    character(*), optional, intent(in)  :: wsynch
!!$    integer(4),   optional, intent(in)  :: PAR_COMM_IN4
!!$    integer(ip)                         :: nsend,nrecv
!!$    integer(ip)                         :: yy_send(2)
!!$    integer(ip)                         :: yy_recv(2)
!!$
!!$    nsend = memor_size(xx_send)
!!$    nrecv = memor_size(xx_recv)
!!$
!!$    if(      nsend == 0 .and. nrecv /= 0 ) then
!!$       call PAR_ALLATOALL_IP(nsend,nrecv,yy_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)
!!$    else if( nsend /= 0 .and. nrecv == 0 ) then
!!$       call PAR_ALLATOALL_IP(nsend,nrecv,xx_send,yy_recv,wherein,wsynch,PAR_COMM_IN4)
!!$    else if( nsend == 0 .and. nrecv == 0 ) then
!!$       return
!!$    else
!!$       call PAR_ALLATOALL_IP(nsend,nrecv,xx_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)
!!$    end if
!!$
!!$  end subroutine PAR_ALLATOALL_IP_1

  subroutine PAR_ALLTOALL_IP_0(nsend,nrecv,xx_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)

    integer(ip),            intent(in)  :: nsend
    integer(ip),            intent(in)  :: nrecv
    integer(ip),            intent(in)  :: xx_send(*)
    integer(ip),            intent(out) :: xx_recv(*)
    character(*), optional, intent(in)  :: wherein
    character(*), optional, intent(in)  :: wsynch
    integer(4),   optional, intent(in)  :: PAR_COMM_IN4
    integer(4)                          :: istat4,nsend4,nrecv4
    integer(4)                          :: PAR_COMM_TO_USE
    logical(lg)                         :: asynch    
    !
    ! Define communicator
    !
    if( present(PAR_COMM_IN4) ) then
       PAR_COMM_TO_USE = PAR_COMM_IN4
    else if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
    end if
    !
    ! BlockinG/non blocking
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    if( asynch ) call runend('ALL TO ALL NOT CODED')
    !
    ! Synchronous Send/receive
    !
    nsend4 = int(nsend,4)
    nrecv4 = int(nrecv,4)
    istat4 = 0

#ifndef MPI_OFF
    call MPI_ALLTOALL(xx_send,nsend4,PAR_INTEGER,xx_recv,nrecv4,PAR_INTEGER,PAR_COMM_TO_USE,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLTOALL_IP_0') 
#endif

  end subroutine PAR_ALLTOALL_IP_0

  subroutine PAR_ALLTOALL_IP_14(xx_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)

    integer(4),             intent(in),    pointer :: xx_send(:)
    integer(4),             intent(inout), pointer :: xx_recv(:)
    character(*), optional, intent(in)             :: wherein
    character(*), optional, intent(in)             :: wsynch
    integer(4),   optional, intent(in)             :: PAR_COMM_IN4
    integer(ip)                                    :: nsend,nrecv,ii
    integer(ip)                                    :: lboun_send,lboun_recv
    integer(4)                                     :: istat4,nsend4,nrecv4
    integer(4)                                     :: PAR_COMM_TO_USE
    integer(4)                                     :: my_rank4,my_size4
    logical(lg)                                    :: asynch

    if( IPARALL ) then
       !
       ! Define communicator
       !
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,my_size4)
       !
       ! BlockinG/non blocking
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if
       if( asynch ) call runend('ALL TO ALL NOT CODED')
       !
       ! Synchronous Send/receive
       !
       nsend      = memory_size(xx_send)/int(my_size4,ip)
       nrecv      = memory_size(xx_recv)/int(my_size4,ip)
       lboun_send = lbound(xx_send,1_ip)
       lboun_recv = lbound(xx_recv,1_ip)
       nsend4     = int(nsend,4)
       nrecv4     = int(nrecv,4)
       istat4     = 0

#ifndef MPI_OFF
       call MPI_ALLTOALL(xx_send(lboun_send:),nsend4,MPI_INTEGER,xx_recv(lboun_recv:),nrecv4,MPI_INTEGER,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLTOALL_IP_14') 
#endif

    else
       if( associated(xx_recv) .and. associated(xx_send) ) then
          do ii = lbound(xx_recv,1),ubound(xx_recv,1)
             xx_recv(ii) = xx_send(ii)
          end do
       end if
    end if

  end subroutine PAR_ALLTOALL_IP_14

  subroutine PAR_ALLTOALL_IP_18(xx_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)

    integer(8),             intent(in),    pointer :: xx_send(:)
    integer(8),             intent(inout), pointer :: xx_recv(:)
    character(*), optional, intent(in)             :: wherein
    character(*), optional, intent(in)             :: wsynch
    integer(4),   optional, intent(in)             :: PAR_COMM_IN4
    integer(ip)                                    :: nsend,nrecv
    integer(ip)                                    :: lboun_send,lboun_recv
    integer(4)                                     :: istat4,nsend4,nrecv4
    integer(4)                                     :: PAR_COMM_TO_USE
    integer(4)                                     :: my_rank4,my_size4
    logical(lg)                                    :: asynch
    !
    ! Define communicator
    !
    if( present(PAR_COMM_IN4) ) then
       PAR_COMM_TO_USE = PAR_COMM_IN4
    else if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,my_size4)
    !
    ! BlockinG/non blocking
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    if( asynch ) call runend('ALL TO ALL NOT CODED')
    !
    ! Synchronous Send/receive
    !
    nsend      = memory_size(xx_send)/int(my_size4,ip)
    nrecv      = memory_size(xx_recv)/int(my_size4,ip)
    lboun_send = lbound(xx_send,1_ip)
    lboun_recv = lbound(xx_recv,1_ip)
    nsend4     = int(nsend,4)
    nrecv4     = int(nrecv,4)
    istat4     = 0

#ifndef MPI_OFF
    call MPI_ALLTOALL(xx_send(lboun_send:),nsend4,MPI_INTEGER8,xx_recv(lboun_recv:),nrecv4,MPI_INTEGER8,PAR_COMM_TO_USE,istat4)
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLTOALL_IP_18') 
#endif

  end subroutine PAR_ALLTOALL_IP_18

  !-----------------------------------------------------------------------
  !>
  !> @author  Ricard Borell
  !> @date    2018-12-04
  !> @brief   All to all v
  !> @details All to all v
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_ALLTOALLV_IP_1(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    integer(ip),  target                        :: send_null(2)
    integer(ip),  target                        :: recv_null(2)
    integer(ip),  pointer                       :: send_tmp(:)
    integer(ip),  pointer                       :: recv_tmp(:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, PAR_INTEGER,  &
            recv_tmp,recvcount4,mpi_rdispls, PAR_INTEGER,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_IP_1

  subroutine PAR_ALLTOALLV_IP_2(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    integer(ip),  target                        :: send_null(2,2)
    integer(ip),  target                        :: recv_null(2,2)
    integer(ip),  pointer                       :: send_tmp(:,:)
    integer(ip),  pointer                       :: recv_tmp(:,:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, PAR_INTEGER,  &
            recv_tmp,recvcount4,mpi_rdispls, PAR_INTEGER,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_IP_2

  subroutine PAR_ALLTOALLV_RP_1(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    real(rp),     target                        :: send_null(2)
    real(rp),     target                        :: recv_null(2)
    real(rp),     pointer                       :: send_tmp(:)
    real(rp),     pointer                       :: recv_tmp(:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, MPI_DOUBLE_PRECISION,  &
            recv_tmp,recvcount4,mpi_rdispls, MPI_DOUBLE_PRECISION,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_RP_1

  subroutine PAR_ALLTOALLV_RP_2(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    real(rp),     target                        :: send_null(2,2)
    real(rp),     target                        :: recv_null(2,2)
    real(rp),     pointer                       :: send_tmp(:,:)
    real(rp),     pointer                       :: recv_tmp(:,:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, MPI_DOUBLE_PRECISION,  &
            recv_tmp,recvcount4,mpi_rdispls, MPI_DOUBLE_PRECISION,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_RP_2

  subroutine PAR_ALLTOALLV_RP_3(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    real(rp),     pointer, intent(in)           :: sendbuf(:,:,:)       !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:,:,:)       !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    real(rp),     target                        :: send_null(2,2,2)
    real(rp),     target                        :: recv_null(2,2,2)
    real(rp),     pointer                       :: send_tmp(:,:,:)
    real(rp),     pointer                       :: recv_tmp(:,:,:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, MPI_DOUBLE_PRECISION,  &
            recv_tmp,recvcount4,mpi_rdispls, MPI_DOUBLE_PRECISION,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_RP_3

  subroutine PAR_ALLTOALLV_RP_4(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    real(rp),     pointer, intent(in)           :: sendbuf(:,:,:,:)       !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:,:,:,:)       !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    real(rp),     target                        :: send_null(2,2,2,2)
    real(rp),     target                        :: recv_null(2,2,2,2)
    real(rp),     pointer                       :: send_tmp(:,:,:,:)
    real(rp),     pointer                       :: recv_tmp(:,:,:,:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, MPI_DOUBLE_PRECISION,  &
            recv_tmp,recvcount4,mpi_rdispls, MPI_DOUBLE_PRECISION,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_RP_4
  !----------------------------------------------------------------------
  !
  ! BROADCAST FOR INTEGERS
  !
  !----------------------------------------------------------------------
  subroutine PAR_BROADCAST_I8_s(xx,wherein,root_rank,PAR_COMM_IN4)
    implicit none
    integer(8),   intent(inout)        :: xx
    character(*),             optional :: wherein
    integer(8),   intent(in), optional :: root_rank
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(8)                         :: n
    integer(4)                         :: root_rank4
    integer(8)                         :: yy(2)
    
    if( ISEQUEN ) return
    if( present(root_rank) ) then
       root_rank4 = int(root_rank,4)
    else
       root_rank4 = 0_4
    end if
    
    yy(1) = xx
    n     = 1
    if( present(PAR_COMM_IN4) ) then
       call PAR_BROADCAST_IP(n,yy,ROOT_RANK4=root_rank4,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_BROADCAST_IP(n,yy,wherein,ROOT_RANK4=root_rank4)
    end if
    xx    = yy(1)
  end subroutine PAR_BROADCAST_I8_s

  subroutine PAR_BROADCAST_I4_s(xx,wherein,PAR_COMM_IN4,root_rank)
    implicit none
    integer(4),   intent(inout)        :: xx
    character(*),             optional :: wherein
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(4),   intent(in), optional :: root_rank
    integer(4)                         :: n
    integer(4)                         :: yy(2)
    if( ISEQUEN ) return
    yy(1) = xx
    n     = 1
    if( present(PAR_COMM_IN4) ) then
       call PAR_BROADCAST_IP(n,yy,ROOT_RANK4=root_rank,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_BROADCAST_IP(n,yy,wherein,ROOT_RANK4=root_rank)
    end if
    xx    = yy(1)
  end subroutine PAR_BROADCAST_I4_s

  subroutine PAR_BROADCAST_I4_1(xx,wherein,PAR_COMM_IN4)

    integer(4),  pointer, intent(inout)         :: xx(:)
    character(*),                      optional :: wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(4)                                  :: n

    if( ISEQUEN ) return

    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if(     present(PAR_COMM_IN4) ) then
             call PAR_BROADCAST_IP(n,xx,PAR_COMM_IN4=PAR_COMM_IN4)
          else if( present(wherein) ) then
             call PAR_BROADCAST_IP(n,xx,wherein)
          else
             call PAR_BROADCAST_I4(n,xx,'IN MY CODE')
          end if
       end if
    end if
  end subroutine PAR_BROADCAST_I4_1

  subroutine PAR_BROADCAST_I8_1(xx,wherein,PAR_COMM_IN4)

    integer(8),   pointer, intent(inout)        :: xx(:)
    character(*),                      optional :: wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(8)                                  :: n
    if( ISEQUEN ) return

    if( associated(xx) ) then
       n = size(xx,KIND=8)
       if( n > 0 ) then
          if(     present(PAR_COMM_IN4) ) then
             call PAR_BROADCAST_IP(n,xx,PAR_COMM_IN4=PAR_COMM_IN4)
          else if( present(wherein) ) then
             call PAR_BROADCAST_IP(n,xx,wherein)
          else
             call PAR_BROADCAST_IP(n,xx,'IN MY CODE')
          end if
       end if
    end if
  end subroutine PAR_BROADCAST_I8_1

  subroutine PAR_BROADCAST_IP_04(n,xx,wherein,PAR_COMM_IN4)

    integer(4),          intent(in)           :: n
    integer(4),          intent(inout)        :: xx(*)
    character(*),                    optional :: wherein
    integer(4),          intent(in), optional :: PAR_COMM_IN4

    if( ISEQUEN ) return

    if( n > 0 ) then
       if(      present(PAR_COMM_IN4) ) then
          call PAR_BROADCAST_IP(n,xx,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( present(wherein) ) then
          call PAR_BROADCAST_IP(n,xx,wherein)
       else
          call PAR_BROADCAST_IP(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_IP_04

  subroutine PAR_BROADCAST_IP_08(n,xx,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),          intent(in)           :: n
    integer(8),          intent(inout)        :: xx(*)
    character(*),                    optional :: wherein
    integer(4),          intent(in), optional :: PAR_COMM_IN4

    if( ISEQUEN .or. n <= 0 ) then
       return
    else
       if(      present(PAR_COMM_IN4) ) then
          call PAR_BROADCAST_IP(n,xx,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( present(wherein) ) then
          call PAR_BROADCAST_IP(n,xx,wherein)
       else
          call PAR_BROADCAST_IP(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_IP_08

  subroutine PAR_BROADCAST_I4(n,xx,wherein,root_rank4,PAR_COMM_IN4)
    implicit none
    integer(4),   intent(in)           :: n
    integer(4),   intent(inout)        :: xx(n)
    character(*),             optional :: wherein
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(4),   intent(in), optional :: root_rank4
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE
    if( ISEQUEN ) return

    if( n > 0 ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(root_rank4) ) then
          call MPI_Bcast(xx(1:n),n4,MPI_INTEGER,root_rank4,PAR_COMM_TO_USE,istat4)
       else
          call MPI_Bcast(xx(1:n),n4,MPI_INTEGER,0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_BROADCAST_I4') 
#endif
    end if
  end subroutine PAR_BROADCAST_I4

  subroutine PAR_BROADCAST_I8(n,xx,wherein,root_rank4,PAR_COMM_IN4)
    implicit none
    integer(8),  intent(in)            :: n
    integer(8),  intent(inout)         :: xx(n)
    character(*),             optional :: wherein
    integer(4),   intent(in), optional :: root_rank4
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE
    if( ISEQUEN ) return

    if( n > 0 ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(root_rank4) ) then
          call MPI_Bcast(xx(1:n),n4,MPI_INTEGER8,root_rank4,PAR_COMM_TO_USE,istat4)
       else
          call MPI_Bcast(xx(1:n),n4,MPI_INTEGER8,0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_BROADCAST_I8') 
#endif
    end if
  end subroutine PAR_BROADCAST_I8


  subroutine PAR_BROADCAST_R8_s(xx,wherein,root_rank)
    implicit none
    real(8),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip),  intent(in), optional :: root_rank
    integer(ip)                 :: n
    real(8)                    :: yy(2)
    yy(1) = xx
    n     = 1
    call PAR_BROADCAST_R8(n,yy,wherein,root_rank)
    xx    = yy(1)
  end subroutine PAR_BROADCAST_R8_s

  subroutine PAR_BROADCAST_R8_1(xx,wherein)
    implicit none
    real(8),     pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: n

    if( ISEQUEN ) return
    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if( present(wherein) ) then
             call PAR_BROADCAST_R8(n,xx,wherein)
          else
             call PAR_BROADCAST_R8(n,xx,'IN MY CODE')
          end if
       end if
    end if

  end subroutine PAR_BROADCAST_R8_1

  subroutine PAR_BROADCAST_R8_2(xx,wherein)
    implicit none
    real(8),     pointer, intent(inout) :: xx(:,:)
    character(*),          optional      :: wherein
    integer(ip)                          :: n

    if( ISEQUEN ) return
    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if( present(wherein) ) then
             call PAR_BROADCAST_R8(n,xx,wherein)
          else
             call PAR_BROADCAST_R8(n,xx,'IN MY CODE')
          end if
       end if
    end if

  end subroutine PAR_BROADCAST_R8_2

  subroutine PAR_BROADCAST_R8_0(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(8),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank

    if( present(wherein) ) then
       if( present(root_rank) ) then
          call PAR_BROADCAST_R8(n,xx,wherein,root_rank)
       else
          call PAR_BROADCAST_R8(n,xx,wherein)
       end if
    else
       if( present(root_rank) ) then
          call PAR_BROADCAST_R8(n,xx,'IN MY CODE',root_rank)
       else
          call PAR_BROADCAST_R8(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_R8_0

  subroutine PAR_BROADCAST_R8_02(n1,n2,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n1
    integer(ip),  intent(in)           :: n2
    real(8),     intent(inout)         :: xx(n1,n2)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank

    if( present(wherein) ) then
       if( present(root_rank) ) then
          call PAR_BROADCAST_R8(n1*n2,xx,wherein,root_rank)
       else
          call PAR_BROADCAST_R8(n1*n2,xx,wherein)
       end if
    else
       if( present(root_rank) ) then
          call PAR_BROADCAST_R8(n1*n2,xx,'IN MY CODE',root_rank)
       else
          call PAR_BROADCAST_R8(n1*n2,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_R8_02

  subroutine PAR_BROADCAST_R8(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(8),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE,root_rank4

    if( ISEQUEN ) return
    if( n > 0 ) then
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
       if( present(root_rank) ) then
          root_rank4 = int(root_rank,4)
       else
          root_rank4 = 0_4
       end if
#ifndef MPI_OFF
       istat4 = 0_4
       call MPI_Bcast(xx(1:n),n4,MPI_DOUBLE_PRECISION,root_rank4,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_BROADCAST_R8')
#endif
    end if
  end subroutine PAR_BROADCAST_R8

  !--4

  subroutine PAR_BROADCAST_R4_s(xx,wherein,root_rank)
    implicit none
    real(4),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip),  intent(in), optional :: root_rank
    integer(ip)                 :: n
    real(4)                    :: yy(2)
    yy(1) = xx
    n     = 1
    call PAR_BROADCAST_R4(n,yy,wherein,root_rank)
    xx    = yy(1)
  end subroutine PAR_BROADCAST_R4_s

  subroutine PAR_BROADCAST_R4_1(xx,wherein)
    implicit none
    real(4),     pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: n

    if( ISEQUEN ) return
    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if( present(wherein) ) then
             call PAR_BROADCAST_R4(n,xx,wherein)
          else
             call PAR_BROADCAST_R4(n,xx,'IN MY CODE')
          end if
       end if
    end if

  end subroutine PAR_BROADCAST_R4_1

  subroutine PAR_BROADCAST_R4_2(xx,wherein)
    implicit none
    real(4),     pointer, intent(inout) :: xx(:,:)
    character(*),          optional      :: wherein
    integer(ip)                          :: n

    if( ISEQUEN ) return
    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if( present(wherein) ) then
             call PAR_BROADCAST_R4(n,xx,wherein)
          else
             call PAR_BROADCAST_R4(n,xx,'IN MY CODE')
          end if
       end if
    end if

  end subroutine PAR_BROADCAST_R4_2

  subroutine PAR_BROADCAST_R4_0(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(4),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank

    if( present(wherein) ) then
       if( present(root_rank) ) then
          call PAR_BROADCAST_R4(n,xx,wherein,root_rank)
       else
          call PAR_BROADCAST_R4(n,xx,wherein)
       end if
    else
       if( present(root_rank) ) then
          call PAR_BROADCAST_R4(n,xx,'IN MY CODE',root_rank)
       else
          call PAR_BROADCAST_R4(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_R4_0

  subroutine PAR_BROADCAST_R4_02(n1,n2,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n1
    integer(ip),  intent(in)           :: n2
    real(4),     intent(inout)         :: xx(n1,n2)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank

    if( present(wherein) ) then
       if( present(root_rank) ) then
          call PAR_BROADCAST_R4(n1*n2,xx,wherein,root_rank)
       else
          call PAR_BROADCAST_R4(n1*n2,xx,wherein)
       end if
    else
       if( present(root_rank) ) then
          call PAR_BROADCAST_R4(n1*n2,xx,'IN MY CODE',root_rank)
       else
          call PAR_BROADCAST_R4(n1*n2,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_R4_02

  subroutine PAR_BROADCAST_R4(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(4),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE,root_rank4

    if( ISEQUEN ) return
    if( n > 0 ) then
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
       if( present(root_rank) ) then
          root_rank4 = int(root_rank,4)
       else
          root_rank4 = 0_4
       end if
#ifndef MPI_OFF
       istat4 = 0_4
       call MPI_Bcast(xx(1:n),n4,MPI_REAL4,root_rank4,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_BROADCAST_R4')
#endif
    end if
  end subroutine PAR_BROADCAST_R4

    !--16

#ifdef __PGI
  subroutine PAR_BROADCAST_R16_s(dummy)
    implicit none
    type(pgi_dummy1), intent(in)                :: dummy
  end subroutine

  subroutine PAR_BROADCAST_R16_1(dummy)
    implicit none
    type(pgi_dummy2), intent(in)                :: dummy
  end subroutine

  subroutine PAR_BROADCAST_R16_0(dummy)
    implicit none
    type(pgi_dummy3), intent(in)                :: dummy
  end subroutine

  subroutine PAR_BROADCAST_R16(dummy)
    implicit none
    type(pgi_dummy5), intent(in)                :: dummy
  end subroutine

#else

  subroutine PAR_BROADCAST_R16_s(xx,wherein)
    implicit none
    real(16),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: n
    real(16)                    :: yy(2)
    yy(1) = xx
    n     = 1
    call PAR_BROADCAST_R16(n,yy,wherein)
    xx    = yy(1)
  end subroutine PAR_BROADCAST_R16_s

  subroutine PAR_BROADCAST_R16_1(xx,wherein)
    implicit none
    real(16),     pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: n

    if( ISEQUEN ) return
    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if( present(wherein) ) then
             call PAR_BROADCAST_R16(n,xx,wherein)
          else
             call PAR_BROADCAST_R16(n,xx,'IN MY CODE')
          end if
       end if
    end if

  end subroutine PAR_BROADCAST_R16_1

  subroutine PAR_BROADCAST_R16_0(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(16),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank

    if( present(wherein) ) then
       if( present(root_rank) ) then
          call PAR_BROADCAST_R16(n,xx,wherein,root_rank)
       else
          call PAR_BROADCAST_R16(n,xx,wherein)
       end if
    else
       if( present(root_rank) ) then
          call PAR_BROADCAST_R16(n,xx,'IN MY CODE',root_rank)
       else
          call PAR_BROADCAST_R16(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_R16_0

  subroutine PAR_BROADCAST_R16(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(16),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE,root_rank4

    if( ISEQUEN ) return
    if( n > 0 ) then
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
       if( present(root_rank) ) then
          root_rank4 = int(root_rank,4)
       else
          root_rank4 = 0_4
       end if
#ifndef MPI_OFF
       istat4 = 0_4
       call runend('PAR_BROADCAST: RP 16 IS NOT COMPATIBLE WITH MPI')
       !call MPI_Bcast(xx(1:n),n4,MPI_REAL4,root_rank4,PAR_COMM_TO_USE,istat4)
       !if( istat4 /= MPI_SUCCESS ) call runend('PAR_BROADCAST_RP: MPI ERROR')
#endif
    end if
  end subroutine PAR_BROADCAST_R16

#endif

  subroutine PAR_BROADCAST_CH(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)           :: n
    character(*), intent(inout)        :: xx
    character(*), intent(in)           :: wherein
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE

    if( IPARALL .and. n > 0 ) then
#ifndef MPI_OFF
       !if( present(wherein) ) then
       !   call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       !else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       !end if
       istat4 = 0_4
       n4 = int(n,4)
       call MPI_Bcast(xx(1:n),n4,MPI_CHARACTER,0_4,PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_BROADCAST_CH')
#endif
    end if

  end subroutine PAR_BROADCAST_CH

  subroutine PAR_BROADCAST_CH_1(xx,wherein)
    implicit none
    character(len=:), intent(inout), pointer :: xx
    character(*),     intent(in)             :: wherein
    integer(4)                               :: istat4,n4,PAR_COMM_TO_USE,n

    if( associated(xx) ) then
       n = len(xx)
       if( IPARALL .and. n > 0 ) then
#ifndef MPI_OFF
          !if( present(wherein) ) then
          !   call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          !else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          !end if
          istat4 = 0_4
          n4 = int(n,4)
          call MPI_Bcast(xx(1:n),n4,MPI_CHARACTER,0_4,PAR_COMM_TO_USE,istat4)
          if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_BROADCAST_CH_1')
#endif
       end if
    end if
    
  end subroutine PAR_BROADCAST_CH_1

  subroutine PAR_BROADCAST_LG_s(xx,wherein)
    implicit none
    logical(lg),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(4)                  :: istat4,PAR_COMM_TO_USE4

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE4)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE4)
       end if
       call MPI_Bcast(xx,1_4,MPI_LOGICAL,0_4,PAR_COMM_TO_USE4,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_BROADCAST_LG_s')
#endif
    end if
  end subroutine PAR_BROADCAST_LG_s

  subroutine PAR_BROADCAST_LG_1(xx,wherein)
    implicit none
    logical(lg),  pointer, intent(inout) :: xx(:)
    character(*), optional               :: wherein
    integer(4)                           :: istat4,PAR_COMM_TO_USE4,n4

    n4 = int(memory_size(xx),4_4)
    if( IPARALL .and. n4 > 0 ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE4)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE4)
       end if
       call MPI_Bcast(xx,n4,MPI_LOGICAL,0_4,PAR_COMM_TO_USE4,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_BROADCAST_LG_1')
#endif
    end if
  end subroutine PAR_BROADCAST_LG_1

  !----------------------------------------------------------------------
  !
  ! PAR_ALL_TO_ALL_ARRAY_OPERATION_IP: FOR INTEGERS
  ! INPUT:  XX(NDOFN) for all slaves
  ! OUTPUT: XX(IDOFN) = XX(IDOFN) if my ranks is the max who have XX(IDOFN) /= 0
  !                   = 0 otherwise
  !         For the master XX(IDOFN) = rank of partition who has XX(IDOFN) /= 0
  !
  !
  !----------------------------------------------------------------------

  subroutine PAR_ALL_TO_ALL_ARRAY_OPERATION_IP(ndofn,xx,what,wherein)
    implicit none
    integer(ip),         intent(in)    :: ndofn
    integer(ip),         intent(inout) :: xx(*)
    character(*),        intent(in)    :: what
    character(*),        optional      :: wherein
    integer(ip)                        :: idofn
    integer(4)                         :: my_rank4
    integer(4)                         :: PAR_COMM_TO_USE4
    integer(ip),         pointer       :: lranks(:)

#ifndef MPI_OFF
    if( IPARALL ) then

       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE4)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE4)
       end if
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE4,my_rank4)

       allocate( lranks(ndofn) )

       if( trim(what) == 'CHOOSE ONLY ONE' ) then
          if( INOTMASTER ) then
             do idofn = 1,ndofn
                if( xx(idofn) > 0 ) then
                   lranks(idofn) = int(my_rank4,ip)
                else
                   lranks(idofn) = 0_ip
                end if
             end do
          end if

          if( present(wherein) ) then
             call PAR_MAX(ndofn,lranks,wherein)
          else
             call PAR_MAX(ndofn,lranks)
          end if

          if( INOTMASTER ) then
             do idofn = 1,ndofn
                if( my_rank4 /= lranks(idofn) ) xx(idofn) = 0
             end do
          else
             do idofn = 1,ndofn
                xx(idofn) = lranks(idofn)
             end do
          end if

       end if

       deallocate( lranks )

    end if
#endif

  end subroutine PAR_ALL_TO_ALL_ARRAY_OPERATION_IP

    !----------------------------------------------------------------------
  !
  ! Gather to rank = 0
  !
  !----------------------------------------------------------------------

  subroutine PAR_GATHER_CHARACTER(character_in,character_out,wherein)
    implicit none
    character(*),           intent(in)    :: character_in
    character(*), pointer,  intent(inout) :: character_out(:)
    character(*), optional, intent(in)    :: wherein
    character(1), pointer                 :: dummi(:)
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: l_character_in
    integer(4)                            :: l_character_out
    integer(4)                            :: istat4

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if

       l_character_in  = len(character_in)

       if( .not. associated(character_out) ) then
          allocate(dummi(l_character_in))
          call MPI_Gather(character_in, l_character_in, MPI_CHARACTER,&
               &          dummi,        l_character_in, MPI_CHARACTER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
          deallocate(dummi)
       else
          l_character_out = len(character_out(1))
          if( l_character_in /= l_character_out ) call runend('PAR_GATHER: WRONG CHARACTER LENGTH')

          call MPI_Gather(character_in, l_character_in, MPI_CHARACTER,&
               &          character_out,l_character_out,MPI_CHARACTER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_CHARACTER')
#endif
    else
       character_out = character_in
    end if

  end subroutine PAR_GATHER_CHARACTER

  subroutine PAR_GATHER_IP_s4(sendbuf,recvbuf,wherein,PAR_COMM_IN4)
    implicit none
    integer(4),             intent(in)    :: sendbuf
    integer(4),   pointer,  intent(inout) :: recvbuf(:)
    character(*), optional, intent(in)    :: wherein
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(4)                            :: dummi(2)
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: istat4,my_rank
    integer(ip)                           :: rl

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,      1_4, MPI_INTEGER,&
               &          dummi,        1_4, MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf,      1_4, MPI_INTEGER,&
               &          recvbuf(rl:), 1_4, MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_IP_s4')
#endif
    else
       recvbuf = sendbuf
    end if

  end subroutine PAR_GATHER_IP_s4

  subroutine PAR_GATHER_IP_s48(sendbuf,recvbuf,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),             intent(in)    :: sendbuf
    integer(4),   pointer,  intent(inout) :: recvbuf(:)
    character(*), optional, intent(in)    :: wherein
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(4)                            :: dummi(2)
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: istat4,my_rank
    integer(ip)                           :: rl
    integer(4)                            :: sendbuf4

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       sendbuf4 = int(sendbuf,4)
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf4,     1_4, MPI_INTEGER,&
               &          dummi,        1_4, MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf4,     1_4, MPI_INTEGER,&
               &          recvbuf(rl:), 1_4, MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_IP_s48')
#endif
    else
       recvbuf = int(sendbuf,4)
    end if

  end subroutine PAR_GATHER_IP_s48

  subroutine PAR_GATHER_IP_s8(sendbuf,recvbuf,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),             intent(in)    :: sendbuf
    integer(8),   pointer,  intent(inout) :: recvbuf(:)
    character(*), optional, intent(in)    :: wherein
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(8)                            :: dummi(2)
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: istat4,my_rank
    integer(ip)                           :: rl

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,      1_4, MPI_INTEGER8,&
               &          dummi,        1_4, MPI_INTEGER8,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf,      1_4, MPI_INTEGER8,&
               &          recvbuf(rl:), 1_4, MPI_INTEGER8,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_IP_s8')
#endif
    else
       recvbuf = sendbuf
    end if

  end subroutine PAR_GATHER_IP_s8

  subroutine PAR_GATHER_IP_14(sendbuf,recvbuf,wherein)
    implicit none
    integer(4),  pointer,  intent(in)    :: sendbuf(:)
    integer(4),  pointer,  intent(inout) :: recvbuf(:)
    character(*),          intent(in)    :: wherein
    integer(4)                           :: dummi(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: sl,rl,lbouns,lbounr

    if( IPARALL ) then
#ifndef MPI_OFF
    istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)

       sendcount4 = memory_size(sendbuf)
       sl = lbound(sendbuf,1)

       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf(sl:), sendcount4,MPI_INTEGER,&
               &          dummi,        sendcount4,MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf(sl:), sendcount4,MPI_INTEGER,&
               &          recvbuf(rl:), sendcount4,MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_IP_14')
#endif
    else
       lbouns  = lbound(sendbuf,1)
       lbounr  = lbound(recvbuf,1)
       sl      = memory_size(sendbuf)
       recvbuf(lbounr:lbounr+sl-1) = sendbuf(lbouns:lbouns+sl-1)
    end if

  end subroutine PAR_GATHER_IP_14

  subroutine PAR_GATHER_IP_18(sendbuf,recvbuf,wherein)
    implicit none
    integer(8),  pointer,  intent(in)    :: sendbuf(:)
    integer(8),  pointer,  intent(inout) :: recvbuf(:)
    character(*),          intent(in)    :: wherein
    integer(8)                           :: dummi(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: sl,rl,lbouns,lbounr

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)

       sendcount4 = memory_size(sendbuf)
       sl = lbound(sendbuf,1)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf(sl:), sendcount4,MPI_INTEGER8,&
               &          dummi,        sendcount4,MPI_INTEGER8,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf(sl:), sendcount4,MPI_INTEGER8,&
               &          recvbuf(rl:), sendcount4,MPI_INTEGER8,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_IP_18')
#endif
    else
       lbouns  = lbound(sendbuf,1)
       lbounr  = lbound(recvbuf,1)
       sl      = memory_size(sendbuf)
       recvbuf(lbounr:lbounr+sl-1) = sendbuf(lbouns:lbouns+sl-1)
    end if

  end subroutine PAR_GATHER_IP_18

  subroutine PAR_GATHER_IP_12(sendbuf,recvbuf,wherein)
    implicit none
    integer(ip),  pointer, intent(in)    :: sendbuf(:)
    integer(ip),  pointer, intent(inout) :: recvbuf(:,:)
    character(*),          intent(in)    :: wherein
    integer(4)                           :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: rl1,rl2,lbouns,sl,ii

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank)

       sendcount4 = memory_size(sendbuf)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,  sendcount4,PAR_INTEGER,&
               &          dummr,    sendcount4,PAR_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl1 = lbound(recvbuf,1)
          rl2 = lbound(recvbuf,2)
          call MPI_Gather(sendbuf,           sendcount4,PAR_INTEGER,&
               &          recvbuf(rl1:,rl2:),sendcount4,PAR_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_IP_12')
#endif
    else
       lbouns = lbound(sendbuf,1)
       rl1    = lbound(recvbuf,1)
       rl2    = lbound(recvbuf,2)
       sl     = memory_size(sendbuf)
       do ii = rl2,ubound(recvbuf,2)
          recvbuf(rl1:rl1+sl-1,ii) = sendbuf(lbouns:lbouns+sl-1)
       end do
    end if

  end subroutine PAR_GATHER_IP_12

  subroutine PAR_GATHER_IP_23(sendbuf,recvbuf,wherein)
    implicit none
    integer(ip),  pointer, intent(in)    :: sendbuf(:,:)
    integer(ip),  pointer, intent(inout) :: recvbuf(:,:,:)
    character(*),          intent(in)    :: wherein
    integer(4)                           :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: rl1,rl2,rl3

#ifndef MPI_OFF
    istat4 = 0_4
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank)

    sendcount4 = memory_size(sendbuf)
    if( .not. associated(recvbuf) ) then
       if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
       call MPI_Gather(sendbuf,  sendcount4,PAR_INTEGER,&
            &          dummr,    sendcount4,PAR_INTEGER,&
            &          0_4,PAR_COMM_TO_USE,istat4)
    else
       rl1 = lbound(recvbuf,1)
       rl2 = lbound(recvbuf,2)
       rl3 = lbound(recvbuf,3)
       call MPI_Gather(sendbuf,                sendcount4,PAR_INTEGER,&
            &          recvbuf(rl1:,rl2:,rl3:),sendcount4,PAR_INTEGER,&
            &          0_4,PAR_COMM_TO_USE,istat4)
    end if
    if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_IP_23')
#endif

  end subroutine PAR_GATHER_IP_23

  subroutine PAR_GATHER_RP_s(sendbuf,recvbuf,wherein)
    implicit none
    real(rp),              intent(in)    :: sendbuf
    real(rp),     pointer, intent(inout) :: recvbuf(:)
    character(*),          intent(in)    :: wherein
    real(rp)                             :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(ip)                          :: rl

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,      1_4, MPI_DOUBLE_PRECISION,&
               &          dummr,        1_4, MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf,      1_4, MPI_DOUBLE_PRECISION,&
               &          recvbuf(rl:), 1_4, MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call runend('PAR_GATHER_IP_s: MPI ERROR')
#endif
    else
       recvbuf = sendbuf
    end if

  end subroutine PAR_GATHER_RP_s

  subroutine PAR_GATHER_RP_1(sendbuf,recvbuf,wherein)
    implicit none
    real(rp),     pointer, intent(in)    :: sendbuf(:)
    real(rp),     pointer, intent(inout) :: recvbuf(:)
    character(*),          intent(in)    :: wherein
    real(rp)                             :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: sl,rl,lbouns,lbounr

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)

       sendcount4 = memory_size(sendbuf)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf, sendcount4,MPI_DOUBLE_PRECISION,&
               &          dummr,   sendcount4,MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf,      sendcount4,MPI_DOUBLE_PRECISION,&
               &          recvbuf(rl:), sendcount4,MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_RP_1')
#endif
    else
       lbouns  = lbound(sendbuf,1)
       lbounr  = lbound(recvbuf,1)
       sl      = memory_size(sendbuf)
       recvbuf(lbounr:lbounr+sl-1) = sendbuf(lbouns:lbouns+sl-1)
    end if

  end subroutine PAR_GATHER_RP_1

  subroutine PAR_GATHER_RP_12(sendbuf,recvbuf,wherein)
    implicit none
    real(rp),     pointer, intent(in)    :: sendbuf(:)
    real(rp),     pointer, intent(inout) :: recvbuf(:,:)
    character(*),          intent(in)    :: wherein
    real(rp)                             :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: rl1,rl2,lbouns,sl

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank)

       sendcount4 = memory_size(sendbuf)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,  sendcount4,MPI_DOUBLE_PRECISION,&
               &          dummr,    sendcount4,MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl1 = lbound(recvbuf,1)
          rl2 = lbound(recvbuf,2)
          call MPI_Gather(sendbuf,           sendcount4,MPI_DOUBLE_PRECISION,&
               &          recvbuf(rl1:,rl2:),sendcount4,MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHER_RP_12')
#endif
    else
       lbouns = lbound(sendbuf,1)
       rl1    = lbound(recvbuf,1)
       rl2    = lbound(recvbuf,2)
       sl     = memory_size(sendbuf)
       recvbuf(rl1:rl1+sl-1,rl2) = sendbuf(lbouns:lbouns+sl-1)
    end if

  end subroutine PAR_GATHER_RP_12

  !-----------------------------------------------------------------------
  !
  !> @brief   Bridge to MPI_ALLGATHERV
  !> @details Bridge to MPI_ALLGATHERV. If the displacement is not
  !>          prescribed the recvbuf are put one after the other
  !>          automatically
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

    subroutine PAR_GATHERV_RP_1(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       sendcount4 = memory_size(sendbuf)
       recvbuf    = sendbuf
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4

       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = PAR_COMM_MY_CODE
       end if
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4 = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_RP_12')
#endif
    end if

  end subroutine PAR_GATHERV_RP_1

  subroutine PAR_GATHERV_RP_0(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),              intent(in)           :: sendbuf(*)           !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in)           :: sendcount4           !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       recvbuf    = sendbuf(1:sendcount4)
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4 = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_RP_0')
#endif
    end if

  end subroutine PAR_GATHERV_RP_0

  subroutine PAR_GATHERV_RP_21(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       sendcount4 = memory_size(sendbuf)
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(kk) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_RP_21')
#endif
    end if

  end subroutine PAR_GATHERV_RP_21

  subroutine PAR_GATHERV_RP_22(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       sendcount4 = memory_size(sendbuf)
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(ii,jj) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4 = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_RP_22')
#endif
    end if

  end subroutine PAR_GATHERV_RP_22

  subroutine PAR_GATHERV_RP_33(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:,:)       !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:,:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk,ll
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       sendcount4 = memory_size(sendbuf)
       kk = lbound(recvbuf,1)
       do ll = lbound(sendbuf,3),ubound(sendbuf,3)
          do jj = lbound(sendbuf,2),ubound(sendbuf,2)
             do ii = lbound(sendbuf,1),ubound(sendbuf,1)
                recvbuf(ii,jj,ll) = sendbuf(ii,jj,ll)
                kk = kk + 1
             end do
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_RP_33')
#endif
    end if

  end subroutine PAR_GATHERV_RP_33

  subroutine PAR_GATHERV_IP_1(sendbuf,recvbuf,recvcount4,wherein,displs4,PAR_COMM_IN4,pp)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4),   pointer, intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(ip), optional, intent(in) :: pp
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if

       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
          sl = lbound(sendbuf,1)
       end if
       if( associated(recvbuf) )    rl  = lbound(recvbuf,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if

          if( present(pp) ) then
             print*,'a=',kfl_paral,associated(sendbuf),size(sendbuf),sendcount4,associated(recvbuf)
          end if
          
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
           end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_IP_1')
    end if
#endif

  end subroutine PAR_GATHERV_IP_1

  subroutine PAR_GATHERV_IP_21(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,rcl,dl
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)
    integer(ip)                                 :: jpart

#ifndef MPI_OFF
    if( IPARALL ) then
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf,1_ip)*memory_size(sendbuf,2_ip)
       end if
       if( associated(recvbuf) )    rl  = lbound(recvbuf,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if
       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_IP_21')
    end if
#endif

  end subroutine PAR_GATHERV_IP_21

  subroutine PAR_GATHERV_IP_22(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       sendcount4 = memory_size(sendbuf)
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(ii,jj) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_IP_22')
#endif
    end if

  end subroutine PAR_GATHERV_IP_22

  subroutine PAR_GATHERV_RP_21_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(kk) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_RP_21_SEND')
#endif
    end if

  end subroutine PAR_GATHERV_RP_21_SEND

  subroutine PAR_GATHERV_RP_22_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(ii,jj) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_RP_22_SEND')
#endif
    end if

  end subroutine PAR_GATHERV_RP_22_SEND

  subroutine PAR_GATHERV_IP_1_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(sendbuf) ) then
          sl = lbound(sendbuf,1)
       end if
       if( associated(recvbuf) )    rl  = lbound(recvbuf,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
           end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_IP_1_SEND')
    end if
#endif

  end subroutine PAR_GATHERV_IP_1_SEND

  subroutine PAR_GATHERV_IP_21_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,rcl,dl
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)
    integer(ip)                                 :: jpart

#ifndef MPI_OFF
    if( IPARALL ) then
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(recvbuf) )    rl  = lbound(recvbuf,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if
       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_IP_21_SEND')
    end if
#endif

  end subroutine PAR_GATHERV_IP_21_SEND


  subroutine PAR_GATHERV_IP_22_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),     pointer, intent(in)        :: sendbuf(:,:)         !< Send buffer
    integer(ip),     pointer, intent(inout)       :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(ii,jj) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_GATHERV_IP_22_SEND')
#endif
    end if

  end subroutine PAR_GATHERV_IP_22_SEND

  !-----------------------------------------------------------------------
  !
  !> @brief   Bridge to MPI_ALLGATHERV
  !> @details Bridge to MPI_ALLGATHERV. If the displacement is not
  !>          prescribed the recvbuf are put one after the other
  !>          automatically
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_ALLGATHERV_IP4(sendbuf,recvbuf,recvcount4,wherein,displs4,PAR_COMM_IN4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart,lbouns,lbounr,lbounc
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( IPARALL ) then
#ifndef MPI_OFF
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHERV_IP4')
#endif
    else
       lbounr = lbound(recvbuf,1)
       lbouns = lbound(sendbuf,1)
       lbounc = lbound(recvcount4,1)
       recvbuf(lbounr:lbounr+recvcount4(lbounc)-1) = sendbuf(lbouns:lbouns+recvcount4(lbounc)-1)
    end if

  end subroutine PAR_ALLGATHERV_IP4

  subroutine PAR_ALLGATHERV_IP8(sendbuf,recvbuf,recvcount8,wherein,displs8,PAR_COMM_IN4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(8),   pointer, intent(in)           :: recvcount8(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(8),   pointer, intent(in), optional :: displs8(:)           !< Displacement
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4),   pointer                       :: displs4(:)
    integer(4),   pointer                       :: recvcount4(:)
    integer(ip)                                 :: jpart
    integer(4)                                  :: lbounr,lbouns,lbounc
    integer(4)                                  :: dumm4

    if( IPARALL ) then
#ifndef MPI_OFF
       sendcount4 = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,4)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       allocate( recvcount4(max(1_ip,size(recvcount8,KIND=ip))) )
       recvcount4 = int(recvcount8,4)
       if( present(displs8) ) then
          if( .not. associated(displs8) ) then
             allocate( displs4(1) )
             displs4 = 0
          else
             allocate( displs4(size(displs8)) )
             displs4 = int(displs8,4)
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,    PAR_INTEGER,&
                  &              recvbuf,recvcount4,displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,        PAR_INTEGER,&
                  &              recvbuf,recvcount4,displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate(displs4)
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount8,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + int(recvcount8(jpart),4)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,       PAR_INTEGER,&
                  &              recvbuf,recvcount4,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,           PAR_INTEGER,&
                  &              recvbuf,recvcount4,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       deallocate( recvcount4 )
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHERV_IP8')
#endif
    else
       lbounr = lbound(recvbuf,1)
       lbouns = lbound(sendbuf,1)
       lbounc = lbound(recvcount8,1)
       dumm4  = int(recvcount8(lbounc),4)
       recvbuf(lbounr:lbounr+dumm4-1) = sendbuf(lbouns:lbouns+dumm4-1)
    end if

  end subroutine PAR_ALLGATHERV_IP8

  subroutine PAR_ALLGATHERV_IP4_2(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHERV_IP4_2')
    end if
#endif

  end subroutine PAR_ALLGATHERV_IP4_2


 subroutine PAR_ALLGATHERV_RP_1(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHERV_RP_1')
    end if
#endif

  end subroutine PAR_ALLGATHERV_RP_1

 subroutine PAR_ALLGATHERV_RP_18(sendbuf,recvbuf,recvcount8,wherein,displs8)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(8),   pointer, intent(in)           :: recvcount8(:)        !< Recv counts
    integer(4),   pointer                       :: recvcount4(:)        !< Recv counts integer 4
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(8),   pointer, intent(in), optional :: displs8(:)           !< Displacement
    integer(4),   pointer                       :: displs4(:)           !< Displacement integer 4
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       allocate( recvcount4(max(1_ip,size(recvcount8,KIND=ip))) )
       recvcount4 = int(recvcount8,4)
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if
      if( present(displs8) ) then
          if( .not. associated(displs8) ) then
             allocate( displs4(1) )
             displs4 = 0
          else
             allocate( displs4(size(displs8)) )
             displs4 = int(displs8,4)
          end if
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHERV_RP_18')
    end if
#endif
  end subroutine PAR_ALLGATHERV_RP_18

  subroutine PAR_ALLGATHERV_RP_2(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHERV_RP_2')
    end if
#endif

  end subroutine PAR_ALLGATHERV_RP_2

  subroutine PAR_ALLGATHERV_RP_3(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:,:)       !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:,:,:)       !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHERV_RP_3')
    end if
#endif

  end subroutine PAR_ALLGATHERV_RP_3

  !-----------------------------------------------------------------------
  !
  !> @brief   Bridge to MPI_ALLGATHER
  !> @details Bridge to MPI_ALLGATHER
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_ALLGATHER_s4(sendbuf,recvbuf,recvcount4_opt,wherein,PAR_COMM_IN4)
    implicit none
    integer(4),            intent(in)           :: sendbuf              !< Send buffer
    integer(4),  pointer,  intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in), optional :: recvcount4_opt       !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Wherein
    integer(4)                                  :: istat4,lbounr
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4)                                  :: recvcount4

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(recvcount4_opt) ) then
          recvcount4 = recvcount4_opt
       else
          recvcount4 = 1_4
       end if
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
          end if
       end if
       sendcount4 = 1_4
       call MPI_ALLGATHER(sendbuf,sendcount4,MPI_INTEGER,&
            &             recvbuf,recvcount4,MPI_INTEGER,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHER_s4')
#endif
    else
       if( associated(recvbuf) ) then
          lbounr = lbound(recvbuf,1)
          recvbuf(lbounr) = sendbuf
       end if
    end if

  end subroutine PAR_ALLGATHER_s4

  subroutine PAR_ALLGATHER_s8(sendbuf,recvbuf,recvcount4_opt,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),            intent(in)           :: sendbuf              !< Send buffer
    integer(8),  pointer,  intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in), optional :: recvcount4_opt       !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Wherein
    integer(4)                                  :: istat4,lbounr
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4)                                  :: recvcount4

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(recvcount4_opt) ) then
          recvcount4 = recvcount4_opt
       else
          recvcount4 = 1_4
       end if
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
          end if
       end if
       sendcount4 = 1_4
       call MPI_ALLGATHER(sendbuf,sendcount4,MPI_INTEGER8,&
            &             recvbuf,recvcount4,MPI_INTEGER8,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHER_s8')
#endif
    else
       if( associated(recvbuf) ) then
          lbounr = lbound(recvbuf,1)
          recvbuf(lbounr) = sendbuf
       end if
    end if

  end subroutine PAR_ALLGATHER_s8

  subroutine PAR_ALLGATHER_IP_14(sendbuf,recvbuf,recvcount4,wherein,PAR_COMM_IN4)
    implicit none
    integer(4),  pointer,  intent(in)             :: sendbuf(:)           !< Send buffer
    integer(4),  pointer,  intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in)             :: recvcount4           !< Recv counts
    character(*),          intent(in),   optional :: wherein              !< Wherein
    integer(4),            intent(in),   optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                    :: istat4
    integer(4)                                    :: lbouns,lbounr
    integer(4)                                    :: PAR_COMM_TO_USE
    integer(4)                                    :: sendcount4
    integer(ip)                                   :: sendnul(2)

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
          end if
       end if
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
          call MPI_ALLGATHER(sendbuf,sendcount4,MPI_INTEGER,&
               &             recvbuf,recvcount4,MPI_INTEGER,&
               &             PAR_COMM_TO_USE,istat4)
       else
          sendcount4 = 0_4
          sendnul    = 0_ip
          call MPI_ALLGATHER(sendnul,sendcount4,MPI_INTEGER,&
               &             recvbuf,recvcount4,MPI_INTEGER,&
               &             PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHER_IP_14')
#endif
    else
       lbounr = lbound(recvbuf,1)
       lbouns = lbound(sendbuf,1)
       recvbuf(lbounr:lbounr+recvcount4-1) = sendbuf(lbouns:lbouns+recvcount4-1)
    end if

  end subroutine PAR_ALLGATHER_IP_14

  subroutine PAR_ALLGATHER_IP_18(sendbuf,recvbuf,recvcount4,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),  pointer,  intent(in)             :: sendbuf(:)           !< Send buffer
    integer(8),  pointer,  intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in)             :: recvcount4           !< Recv counts
    character(*),          intent(in),   optional :: wherein              !< Wherein
    integer(4),            intent(in),   optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                    :: istat4
    integer(4)                                    :: lbouns,lbounr
    integer(4)                                    :: PAR_COMM_TO_USE
    integer(4)                                    :: sendcount4
    integer(ip)                                   :: sendnul(2)

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
          end if
       end if
       if( associated(sendbuf) ) then
          sendcount4 = memory_size(sendbuf)
          call MPI_ALLGATHER(sendbuf,sendcount4,MPI_INTEGER8,&
               &             recvbuf,recvcount4,MPI_INTEGER8,&
               &             PAR_COMM_TO_USE,istat4)
       else
          sendcount4 = 0_4
          sendnul    = 0_ip
          call MPI_ALLGATHER(sendnul,sendcount4,MPI_INTEGER8,&
               &             recvbuf,recvcount4,MPI_INTEGER8,&
               &             PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHER_IP_18')
#endif
    else
       lbounr = lbound(recvbuf,1)
       lbouns = lbound(sendbuf,1)
       recvbuf(lbounr:lbounr+recvcount4-1) = sendbuf(lbouns:lbouns+recvcount4-1)
    end if

  end subroutine PAR_ALLGATHER_IP_18
  
  subroutine PAR_ALLGATHER_RP(sendbuf,recvbuf,wherein)
    implicit none
    real(rp),              intent(in)  :: sendbuf             !< Send buffer
    real(rp),              intent(out) :: recvbuf(*)           !< Recv buffer
    character(*),          intent(in)  :: wherein                !< Wherein
    integer(4)                         :: istat4
    integer(4)                         :: PAR_COMM_TO_USE

#ifndef MPI_OFF
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_ALLGATHER(sendbuf,1_4,MPI_DOUBLE_PRECISION,&
            &             recvbuf,1_4,MPI_DOUBLE_PRECISION,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHER_IP_14')
    end if
#endif

  end subroutine PAR_ALLGATHER_RP

  subroutine PAR_ALLGATHER_RP_0(sendbuf,recvbuf,recvcount4,wherein)
    implicit none
    real(rp),              intent(in)  :: sendbuf(*)           !< Send buffer
    real(rp),              intent(out) :: recvbuf(*)           !< Recv buffer
    integer(4),            intent(in)  :: recvcount4           !< Recv counts
    character(*),          intent(in)  :: wherein                !< Wherein
    integer(4)                         :: istat4
    integer(4)                         :: PAR_COMM_TO_USE

#ifndef MPI_OFF
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_ALLGATHER(sendbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             recvbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHER_RP_0')
    end if
#endif

  end subroutine PAR_ALLGATHER_RP_0

  subroutine PAR_ALLGATHER_RP_2(sendbuf,recvbuf,recvcount4,wherein)
    implicit none
    real(rp),     pointer, intent(in)    :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout) :: recvbuf(:,:)         !< Recv buffer
    integer(4),            intent(in)    :: recvcount4           !< Recv counts
    character(*),          intent(in)    :: wherein                !< Wherein
    integer(4)                           :: istat4
    integer(4)                           :: PAR_COMM_TO_USE

#ifndef MPI_OFF
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_ALLGATHER(sendbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             recvbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHER_RP_2')
    end if
#endif

  end subroutine PAR_ALLGATHER_RP_2

  subroutine PAR_ALLGATHER_RP_02(sendbuf,recvbuf,recvcount4,wherein,PAR_COMM_IN4)
    implicit none
    real(rp),              intent(in)    :: sendbuf(*)           !< Send buffer
    real(rp),     pointer, intent(inout) :: recvbuf(:,:)         !< Recv buffer
    integer(4),            intent(in)    :: recvcount4           !< Recv counts
    character(*), optional,intent(in)    :: wherein              !< Wherein
    integer(4),   optional,intent(in)    :: PAR_COMM_IN4         !< Communicator
    integer(4)                           :: istat4,lbounr,lbouns
    integer(4)                           :: PAR_COMM_TO_USE

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       end if
       call MPI_ALLGATHER(sendbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             recvbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHER_RP_02')
#endif
    else
       lbounr = lbound(recvbuf,2)
       recvbuf(1:recvcount4,lbounr) = sendbuf(1:recvcount4)
    end if

  end subroutine PAR_ALLGATHER_RP_02

  subroutine PAR_ALLGATHER_LG(sendbuf,recvbuf,recvcount4,wherein)
    implicit none
    logical(lg),           intent(in)  :: sendbuf              !< Send buffer
    logical(lg),  pointer, intent(inout) :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in)  :: recvcount4           !< Recv counts
    character(*),          intent(in)  :: wherein                !< Wherein
    integer(4)                         :: istat4
    integer(4)                         :: PAR_COMM_TO_USE
    integer(4)                         :: sendcount4

    if( IPARALL ) then
#ifndef MPI_OFF
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       sendcount4 = 1_4
       call MPI_ALLGATHER(sendbuf,sendcount4,MPI_LOGICAL,&
            &             recvbuf,recvcount4,MPI_LOGICAL,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_ALLGATHER_LG')
#endif
    end if

  end subroutine PAR_ALLGATHER_LG

  !----------------------------------------------------------------------
  !
  ! Scatter from rank = 0
  !
  !----------------------------------------------------------------------

  subroutine PAR_SCATTER_IP_s(xx_in,xx_out,wherein)
    implicit none
    integer(ip),  pointer, intent(in)  :: xx_in(:)
    integer(ip),           intent(out) :: xx_out
    character(*),          intent(in)  :: wherein
    integer(ip)                        :: dummi(2)
    integer(4)                         :: PAR_COMM_TO_USE
    integer(4)                         :: istat4

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( .not. associated(xx_in) ) then
          call MPI_Scatter(dummi,1_4,PAR_INTEGER,xx_out,1_4,PAR_INTEGER,0_4,PAR_COMM_TO_USE,istat4)
       else
          call MPI_Scatter(xx_in,1_4,PAR_INTEGER,xx_out,1_4,PAR_INTEGER,0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTER_IP_s')
    end if
#endif

  end subroutine PAR_SCATTER_IP_s

  subroutine PAR_SCATTERV_IP_1(sendbuf,recvbuf,sendcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4
       recvcount4      = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          recvcount4 = memory_size(recvbuf)
          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

!       if( associated(sendcount4) ) then
       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,              PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                      PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,             PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                     PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTER_IP_1')

    end if

#endif

  end subroutine PAR_SCATTERV_IP_1

  subroutine PAR_SCATTERV_IP_2(sendbuf,recvbuf,sendcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)           !< Send buffer
    integer(ip),  pointer, intent(inout)          :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4
       recvcount4      = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          recvcount4 = memory_size(recvbuf)
          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

!       if( associated(sendcount4) ) then
       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,              PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                      PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,             PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                     PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTER_IP_2')

    end if

#endif

  end subroutine PAR_SCATTERV_IP_2

  subroutine PAR_SCATTERV_RP_0(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),              intent(inout)        :: recvbuf(*)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),            intent(in)           :: recvcount4           !< Recv count
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),   pointer                       :: my_displs4(:)
    real(rp)                                    :: sendbuf_tmp(2)
    real(rp)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       rl = 1
       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

!       if( associated(sendcount4) ) then
       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTER_RP_0')

    end if

#endif

  end subroutine PAR_SCATTERV_RP_0

  subroutine PAR_SCATTERV_RP_1(sendbuf,recvbuf,sendcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(rp)                                    :: sendbuf_tmp(2)
    real(rp)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4
       recvcount4      = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          recvcount4 = memory_size(recvbuf)
          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

!       if( associated(sendcount4) ) then
       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTER_RP_1')

    end if

#endif

  end subroutine PAR_SCATTERV_RP_1

  subroutine PAR_SCATTERV_RP_2(sendbuf,recvbuf,sendcount4,wherein,displs4)
    implicit none
    real(rp),  pointer, intent(in)              :: sendbuf(:,:)           !< Send buffer
    real(rp),  pointer, intent(inout)             :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)          !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)             !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(rp)                                    :: sendbuf_tmp(2)
    real(rp)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),    pointer                      :: sendcount4_tmp(:)
    integer(4),    target                       :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4
       recvcount4      = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          recvcount4 = memory_size(recvbuf)
          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTER_RP_2')

    end if

#endif

  end subroutine PAR_SCATTERV_RP_2

  subroutine PAR_SCATTERV_IP_1_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    integer(4),            intent(in)           :: recvcount4
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,              PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                      PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,             PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                     PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTERV_IP_1_RCV')

    end if

#endif

  end subroutine PAR_SCATTERV_IP_1_RCV

  subroutine PAR_SCATTERV_IP_2_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),            intent(in)           :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,              PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                      PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,             PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                     PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTERV_IP_2_RCV')

    end if

#endif

  end subroutine PAR_SCATTERV_IP_2_RCV

  subroutine PAR_SCATTERV_RP_1_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),  pointer, intent(in)              :: sendbuf(:)           !< Send buffer
    real(rp),  pointer, intent(inout)           :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),    intent(in)                   :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(rp)                                    :: sendbuf_tmp(2)
    real(rp)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),    pointer                      :: sendcount4_tmp(:)
    integer(4),    target                       :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTERV_RP_1_RCV')

    end if

#endif

  end subroutine PAR_SCATTERV_RP_1_RCV

  subroutine PAR_SCATTERV_R8_2_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(8),  pointer, intent(in)               :: sendbuf(:,:)           !< Send buffer
    real(8),  pointer, intent(inout)            :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)          !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)             !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),     intent(in)                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(8)                                    :: sendbuf_tmp(2)
    real(8)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),    pointer                      :: sendcount4_tmp(:)
    integer(4),    target                       :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTERV_R8_2_RCV')

    end if

#endif

  end subroutine PAR_SCATTERV_R8_2_RCV

  subroutine PAR_SCATTERV_R4_2_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(4),  pointer, intent(in)               :: sendbuf(:,:)           !< Send buffer
    real(4),  pointer, intent(inout)            :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)          !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)             !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),     intent(in)                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(4)                                    :: sendbuf_tmp(2)
    real(4)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),    pointer                      :: sendcount4_tmp(:)
    integer(4),    target                       :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_REAL4,&
                     &            recvbuf_tmp, recvcount4,              MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_REAL4,&
                     &            recvbuf_tmp, recvcount4,                  MPI_REAL4,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_REAL4,&
                     &            recvbuf, recvcount4,                  MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_REAL4,&
                     &            recvbuf, recvcount4,                      MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_REAL4,&
                     &            recvbuf_tmp, recvcount4,             MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_REAL4,&
                     &            recvbuf_tmp, recvcount4,                 MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_REAL4,&
                     &            recvbuf, recvcount4,                 MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_REAL4,&
                     &            recvbuf, recvcount4,                     MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= MPI_SUCCESS ) call PAR_MPI_RUNEND(istat4,'PAR_SCATTERV_R4_2_RCV')

    end if

#endif

  end subroutine PAR_SCATTERV_R4_2_RCV

#ifdef __PGI
  subroutine PAR_SCATTERV_R16_2_RCV(dummy)
    implicit none
    type(pgi_dummy1), intent(in)                :: dummy
  end subroutine PAR_SCATTERV_R16_2_RCV
#else

  subroutine PAR_SCATTERV_R16_2_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(16),  pointer, intent(in)              :: sendbuf(:,:)           !< Send buffer
    real(16),  pointer, intent(inout)           :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)          !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)             !< Displacement
    integer(4),     intent(in)                  :: recvcount4

#ifndef MPI_OFF
    if( IPARALL ) then
        call runend('PAR_BROADCAST: RP 16 IS NOT COMPATIBLE WITH MPI')
    end if

#endif

  end subroutine PAR_SCATTERV_R16_2_RCV

#endif


  !----------------------------------------------------------------------
  !
  ! WAITALL
  !
  !----------------------------------------------------------------------

  subroutine PAR_WAITALL()
    integer(4) :: istat4
    integer(4) :: count4
#ifndef MPI_OFF
    integer(4) :: status4(MPI_STATUS_SIZE,1_ip)
#endif

#ifndef MPI_OFF
    if( IPARALL ) then
       count4 = 1_4
       call MPI_WAITALL(count4,ireq41,status4,istat4)
       deallocate(xx_non_blocking)
    end if
#endif
  end subroutine PAR_WAITALL

end module mod_communications_global
!> @}
