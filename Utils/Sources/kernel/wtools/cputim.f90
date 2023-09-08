!-----------------------------------------------------------------------
!> @addtogroup CPU_Time
!> @{
!> @file    cputim.f90
!> @author  houzeaux
!> @date    2018-12-30
!> @brief   Compute current time
!> @details This routine finds out the CPU time in seconds
!-----------------------------------------------------------------------

subroutine cputim(rtime)

  use def_kintyp, only  :  rp
  use def_master, only  :  rate_time
  implicit none
  real(rp), intent(out) :: rtime
  real(8)               :: rtim8
  real(4)               :: rtim4,elap4(2),etime
  integer(4)            :: itim4
  integer(8)            :: itim8

  if( 1 == 1 ) then
    
     call system_clock(itim8)
     rtime = real(itim8,rp) * rate_time

  else if( 1 == 2 ) then

     call cpu_time(rtime)

  else if( 1 == 3 ) then
     !
     ! This method requires linking with -lrt
     ! Commented temporarly (discussion is required)
     ! ! Using C Wall time routine
     ! call wall_time(rtim8)
     ! ! Just in case rp is not 8
     ! rtime = rtim8
  end if

end subroutine cputim
