subroutine cputim(rtime)
  !-----------------------------------------------------------------------
  ! NAME
  !    cputim
  ! DESCRIPTION
  !    This routine finds out the CPU time in seconds
  !-----------------------------------------------------------------------
  implicit none
  real    , intent(out) :: rtime
  real(4)               :: rtim4,elap4(2),etime
  real                  :: OMP_GET_WTIME

!  rtim4 = etime(elap4)
!  rtime = real(rtim4)
  call cpu_time(rtime)

!  rtime= OMP_GET_WTIME()

end subroutine cputim
