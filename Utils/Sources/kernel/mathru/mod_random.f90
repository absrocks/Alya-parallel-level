!____________________________________________________________________________
! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!***********************************************************************
! Fortran version rewritten as an F90 module and mt state saving and getting
! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999

module mod_random
  !
  use def_kintyp,         only : rp,ip,lg
#ifndef I_AM_NOT_ALYA  
  use mod_communications, only : PAR_BROADCAST
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
#endif
  implicit none

  private
  !
  integer,      parameter              :: ip4 = 4
  ! Default seed
  integer(ip4), parameter              :: defaultsd = 4357
  ! Period parameters
  integer(ip4), parameter              :: N = 624, N1 = N + 1
  ! the array for the state vector
  integer(ip4), save, dimension(0:N-1) :: mt
  integer(ip4), save                   :: mti = N1

  public :: random_grnd
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

  subroutine random_initialization()

    mti = N1

  end subroutine random_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-08-24
  !> @brief   Seeds
  !> @details setting initial seeds to mt[N] using
  !>          the generator Line 25 of Table 1 in
  !>          [KNUTH 1981, The Art of Computer Programming
  !>          Vol. 2 (2nd Ed.), pp102]
  !>          mt(mti) = iand(69069 * mt(mti-1_ip4),-1_ip4) was changed
  !>          to avoid overflows
  !> 
  !-----------------------------------------------------------------------

  subroutine random_sgrnd(seed)
    implicit none

    integer(ip4), intent(in) :: seed

    mt(0) = iand(seed,-1_ip4)
    do mti = 1,N-1
       mt(mti) = int(iand(int(69069,8) * int(mt(mti-1_ip4),8),-1_8),ip4)
    end do
    
  end subroutine random_sgrnd

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-08-24
  !> @brief   Random number generator 
  !> @details Random number generator
  !>          mag01(x) = x * MATA for x=0,1
  !> 
  !-----------------------------------------------------------------------
  
  real(rp) function random_grnd(UNIQUE_SEED,PAR_COMM_IN4,PAR_RANK4)
    
    logical(lg),  intent(in), optional :: UNIQUE_SEED
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(4),   intent(in), optional :: PAR_RANK4
    ! Period parameters
    integer(ip4), parameter            :: M = 397, MATA  = -1727483681 ! constant vector a
    integer(ip4), parameter            :: LMASK =  2147483647          ! least significant r bits
    integer(ip4), parameter            :: UMASK = -LMASK - 1           ! most significant w-r bits
    ! Tempering parameters
    integer(ip4), parameter            :: TMASKB= -1658038656
    integer(ip4), parameter            :: TMASKC= -272236544
    integer(ip4), save                 :: mag01(0:1) = (/ 0, MATA /)
    integer(ip4)                       :: kk,y
    integer(ip4)                       :: TSHFTU,TSHFTS,TSHFTT,TSHFTL

    TSHFTU(y) = ishft(y,-11_ip4)
    TSHFTS(y) = ishft(y,  7_ip4)
    TSHFTT(y) = ishft(y, 15_ip4)
    TSHFTL(y) = ishft(y,-18_ip4)
    
    if( mti >= N ) then
       
       ! generate N words at one time
       
       if( mti == N+1 ) then
          call random_generate_seed(UNIQUE_SEED,PAR_COMM_IN4,PAR_RANK4)          
       end if

       do kk = 0,N-M-1
          y      = ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk) = ieor(ieor(mt(kk+M),ishft(y,-1_ip4)),mag01(iand(y,1_ip4)))
       enddo
       do kk = N-M,N-2
          y      = ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
          mt(kk) = ieor(ieor(mt(kk+(M-N)),ishft(y,-1_ip4)),mag01(iand(y,1_ip4)))
       enddo
       y       = ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
       mt(N-1) = ieor(ieor(mt(M-1),ishft(y,-1_ip4)),mag01(iand(y,1_ip4)))
       mti     = 0
    endif

    y   = mt(mti)
    mti = mti + 1_ip4
    y   = ieor(y,TSHFTU(y))
    y   = ieor(y,iand(TSHFTS(y),TMASKB))
    y   = ieor(y,iand(TSHFTT(y),TMASKC))
    y   = ieor(y,TSHFTL(y))
    if( y < 0_ip4 )  then
       random_grnd = (real(y,rp)+2.0_rp**32)/(2.0_rp**32-1.0_rp)
    else
       random_grnd = real(y,rp)/(2.0_rp**32-1.0_rp)
    endif

  end function random_grnd

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
    integer(ip4)                       :: ij,nn
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
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_IN4,my_rank4)
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
       call random_seed(size = nn)
       allocate(rand_seed(INT(nn,ip4))) 
       !CALL SYSTEM_CLOCK(COUNT=clock)
       CALL DATE_AND_TIME(VALUES=values) ! Using only miliseconds of the time with values(8)
       rand_seed    = values(8) + 37_ip4 * (/ (ij - 1_ip4, ij = 1, int(nn,ip4)) /)
       !
       ! Include my rank
       !
       !rand_seed(1) = XOR(rand_seed(1),INT(ABS(kfl_paral),ip4))
       rand_seed(1) = IEOR(rand_seed(1),abs(my_rank4+1099279_4))
       !
       !!! OJO pruebo sin kfl_paral  rand_seed(1) = XOR(rand_seed(1),INT(ABS(1099279),ip4))
       ! if sgrnd() has not been called,
       !
       call random_sgrnd( int(rand_seed(1),ip4) )
       !call random_sgrnd( defaultsd  ) a default initial seed is used
    end if
    !
    ! Broacast seed
    !
    if( if_unique_seed .and. my_rank4 >= 0_4 .and. present(PAR_COMM_IN4) ) then
#ifndef I_AM_NOT_ALYA  
       call PAR_BROADCAST(nn,PAR_COMM_IN4=PAR_COMM_IN4)
       call PAR_BROADCAST(N,mt(0:N-1),PAR_COMM_IN4=PAR_COMM_IN4)
#else
       call runend('MOD_RANDOM: YOU SHOULD BE ALYA TO EMPLOY A UNIQUE SEED')
#endif
    end if

    if( allocated(rand_seed) ) deallocate(rand_seed)

  end subroutine random_generate_seed

end module mod_random
