subroutine par_barrie()
!------------------------------------------------------------------------
!****f* Parall/par_send
! NAME
!    par_barrie
! DESCRIPTION
!    MPI Barrier
! OUTPUT
!   
! USED BY
!    Parall
!***
!------------------------------------------------------------------------
  use def_parall
  use def_master
  use mod_iofile
  use mod_parall, only : PAR_COMM_MY_CODE4
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h'
#endif
  integer(4) :: istat

#ifdef MPI_OFF
#else
  call MPI_Barrier( PAR_COMM_MY_CODE4, istat )
#endif

end subroutine par_barrie
