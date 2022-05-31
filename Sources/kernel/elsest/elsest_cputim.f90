subroutine elsest_cputim(rtime)
  !-----------------------------------------------------------------------
  !****f* elsest_cputim
  ! NAME
  !    nsi_elmope
  ! DESCRIPTION
  !    Returns the CPU time in seconds
  !    corrections by Alistair Hart <ahart@cray.com> to avoid non standard etime
  ! OUTPUT
  !    rtime
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
#ifdef _OPENMP
  use omp_lib
#endif
  use def_elsest, only   : rp

!!$ CRAY change
!!$   etime() is not part of the Fortran standard
!!$   Some compilers support it, but it would be better to use
!!$   something more standard. 
!!$   In the meantime, we will use MPI_WTIME()
!!$   Macro _CRAYFTN is automatically defined if we are using
!!$   the Cray compiler
#ifdef _CRAYFTN
  use MPI
#endif
!!$ END CRAY change

  implicit none
  real(rp), intent(out) :: rtime
  real(4)               :: elap4(2),etime

#ifdef _OPENMP
  rtime = OMP_GET_WTIME()
#else
!!$ CRAY change
#ifdef _CRAYFTN
  rtime = MPI_WTIME()
#else
  call cpu_time(rtime)
#endif
!!$ END CRAY change
#endif

end subroutine elsest_cputim
