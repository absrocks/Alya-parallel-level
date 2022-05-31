!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    mod_random.f90
!> @author  houzeaux
!> @date    2020-05-11
!> @brief   Random number generator
!> @details Random number generator
!>         
!-----------------------------------------------------------------------

module mod_random
 
  use def_kintyp_basic,   only : rp,ip,lg
#ifndef I_AM_NOT_ALYA  
  use mod_communications, only : PAR_BROADCAST
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use def_master,         only : ISEQUEN 
#endif
  implicit none

  private
 
  integer,      parameter  :: ip4 = 4
  integer(ip4), parameter  :: defaultsd = 4357
  integer(ip4), parameter  :: nn = 624, n1 = nn + 1
  integer(ip4)             :: mt(0:nn-1)
  integer(ip4)             :: mti = n1
  integer(ip4)             :: kfl_different_seed = 1

  public :: random_generate_number
  public :: random_initialization

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-08-24
  !> @brief   Initialization 
  !> @details Initialization 
  !> 
  !-----------------------------------------------------------------------

  subroutine random_initialization(kfl_randseed)
    integer(ip),  intent(in), optional :: kfl_randseed

    if (present(kfl_randseed)) then
        if (kfl_randseed == 0) then
            kfl_different_seed = 0
        endif
    endif

    mti = n1

  end subroutine random_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-08-24
  !> @brief   Seeds
  !> @details setting initial seeds to mt[N]
  !> 
  !-----------------------------------------------------------------------

  subroutine random_seed_table(seed)
    implicit none

    integer(ip4), intent(in) :: seed

    mt(0) = iand(seed,-1_ip4)
    do mti = 1,nn-1
       mt(mti) = int(iand(int(69069,8) * int(mt(mti-1_ip4),8),-1_8),ip4)
    end do
    
  end subroutine random_seed_table

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-08-24
  !> @brief   Random number generator 
  !> @details Random number generator
  !> 
  !-----------------------------------------------------------------------  
  
  real(rp) function random_generate_number(UNIQUE_SEED,PAR_COMM_IN4,PAR_RANK4)
    
    logical(lg),  intent(in), optional :: UNIQUE_SEED
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(4),   intent(in), optional :: PAR_RANK4
    integer(ip4), parameter            :: M = 397, MATA  = -1727483681 ! constant vector a
    integer(ip4), parameter            :: maskl =  2147483647          ! least significant r bits
    integer(ip4), parameter            :: masku = -maskl - 1           ! most significant w-r bits
    integer(ip4), parameter            :: maskb= -1658038656
    integer(ip4), parameter            :: maskc= -272236544
    integer(ip4), save                 :: mag01(0:1) = (/ 0, MATA /)
    integer(ip4)                       :: kk,y
    integer(ip4)                       :: tshftu,tshfts,tshftt,tshftl

    tshftu(y) = ishft(y,-11_ip4)
    tshfts(y) = ishft(y,  7_ip4)
    tshftt(y) = ishft(y, 15_ip4)
    tshftl(y) = ishft(y,-18_ip4)
    
    if( mti >= nn ) then
       
       ! generate nn words at one time
       
       if( mti == nn+1 ) then
          call random_generate_seed(UNIQUE_SEED,PAR_COMM_IN4,PAR_RANK4)          
       end if

       do kk = 0,nn-M-1
          y      = ior(iand(mt(kk),masku),iand(mt(kk+1),maskl))
          mt(kk) = ieor(ieor(mt(kk+M),ishft(y,-1_ip4)),mag01(iand(y,1_ip4)))
       enddo
       do kk = nn-M,nn-2
          y      = ior(iand(mt(kk),masku),iand(mt(kk+1),maskl))
          mt(kk) = ieor(ieor(mt(kk+(M-nn)),ishft(y,-1_ip4)),mag01(iand(y,1_ip4)))
       enddo
       y       = ior(iand(mt(nn-1),masku),iand(mt(0),maskl))
       mt(nn-1) = ieor(ieor(mt(M-1),ishft(y,-1_ip4)),mag01(iand(y,1_ip4)))
       mti     = 0
    endif

    y   = mt(mti)
    mti = mti + 1_ip4
    y   = ieor(y,tshftu(y))
    y   = ieor(y,iand(tshfts(y),maskb))
    y   = ieor(y,iand(tshftt(y),maskc))
    y   = ieor(y,tshftl(y))
    if( y < 0_ip4 )  then
       random_generate_number = (real(y,rp)+2.0_rp**32)/(2.0_rp**32-1.0_rp)
    else
       random_generate_number = real(y,rp)/(2.0_rp**32-1.0_rp)
    endif

  end function random_generate_number

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-08-23
  !> @brief   Seed generator
  !> @details When a unique seed is used, only MPI rank=0 generates the
  !>          seed and broadcast it to antoher ranks.
  !>          This is used when all ranks shoulds generate the same
  !>          random numbers
  !> 
  !-----------------------------------------------------------------------

  subroutine random_generate_seed(UNIQUE_SEED,PAR_COMM_IN4,PAR_RANK4)

    logical(lg),  intent(in), optional :: UNIQUE_SEED   !< If a unique seed should be used in the communicator
    integer(4),   intent(in), optional :: PAR_COMM_IN4  !< MPI communicator
    integer(4),   intent(in), optional :: PAR_RANK4     !< MPI rank
    integer(4)                         :: my_rank4       
    logical(lg)                        :: if_unique_seed
    integer(ip4), allocatable          :: rand_seed(:)
    integer(ip4)                       :: ij,nsize
    integer(ip4), dimension(8)         :: values
    !
    ! If a unique seed should be shared in the communicator
    !
    if_unique_seed = .false.
    if( present(UNIQUE_SEED) ) if_unique_seed = UNIQUE_SEED
    !
    ! Are we in sequential MY_RANK4=-1/0 or parallel mode MY_RANK4 >=0 ?
    !
    my_rank4 = -1_4    
    if(      present(PAR_RANK4)    ) then
       my_rank4 = PAR_RANK4
    else if( present(PAR_COMM_IN4) ) then
#ifndef I_AM_NOT_ALYA  
       if (ISEQUEN) then
           my_rank4 = -1_4
       else
           call PAR_COMM_RANK_AND_SIZE(PAR_COMM_IN4,my_rank4)
       endif
#else
       call runend('MOD_RANDOM: YOU SHOULD BE ALYA TO EMPLOY THIS OPTION... JUST SEND THE MPI RANK ALSO')
#endif
    end if
    !
    ! Generate seed
    !
    if( my_rank4 <= 0_4 .or. ( .not. if_unique_seed ) ) then
       !
       ! generate a random seed using date and time
       !
       call random_seed(size = nsize)
       allocate(rand_seed(INT(nsize,ip4))) 

       if (kfl_different_seed /= 0)  then
          !CALL SYSTEM_CLOCK(COUNT=clock)
          CALL DATE_AND_TIME(VALUES=values) ! Using only miliseconds of the time with values(8)
          rand_seed    = values(8) + 37_ip4 * (/ (ij - 1_ip4, ij = 1, int(nsize,ip4)) /)
          !
          ! Include my rank
          !
          !rand_seed(1) = XOR(rand_seed(1),INT(ABS(kfl_paral),ip4))
          rand_seed(1) = IEOR(rand_seed(1),abs(my_rank4+1099279_4))
       else
          rand_seed    = (/ (1_ip4, ij = 1, int(nsize,ip4)) /)
       endif
       !
       !!! OJO pruebo sin kfl_paral  rand_seed(1) = XOR(rand_seed(1),INT(ABS(1099279),ip4))
       ! if sgrnd() has not been called,
       !
       call random_seed_table( int(rand_seed(1),ip4) )
       !call random_seed_table( defaultsd  ) a default initial seed is used
    end if
    !
    ! Broacast seed
    !
    if( if_unique_seed .and. my_rank4 >= 0_4 .and. present(PAR_COMM_IN4) ) then

#ifndef I_AM_NOT_ALYA  
       call PAR_BROADCAST(nsize,PAR_COMM_IN4=PAR_COMM_IN4)
       call PAR_BROADCAST(nn,mt(0:nn-1),PAR_COMM_IN4=PAR_COMM_IN4)
#else
       call runend('MOD_RANDOM: YOU SHOULD BE ALYA TO EMPLOY A UNIQUE SEED')
#endif

    end if

    if( allocated(rand_seed) ) deallocate(rand_seed)

  end subroutine random_generate_seed

end module mod_random
