!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_alya2signal.f90
!> @author  houzeaux
!> @date    2020-02-25
!> @brief   Signal handling
!> @details Module to treat signal handling... which is not part of
!>          the standard. See https://en.wikipedia.org/wiki/Signal_(IPC)
!-----------------------------------------------------------------------

module mod_alya2signal

  use def_master, only : ITASK_WRITE_RESTART
  use def_master, only : kfl_preli
  use def_master, only : kfl_timei
  use def_master, only : kfl_stop
  implicit none
  
  integer(4), parameter :: SIGINT  =  2_4 ! Interrupt from keyboard Ctrl-C 
  integer(4), parameter :: SIGTERM = 15_4 ! Termination signal             
  integer(4), parameter :: SIGUSR1 = 10_4 ! SIGUSR1
  integer(4), parameter :: SIGUSR2 = 12_4 ! SIGUSR2
  
contains

  subroutine alya2signal

    
    integer(4) :: status
    integer(4) :: flag

    flag=-1_4

#ifdef ALYA_SIGNAL
#if defined __INTEL_COMPILER || __PGI
    status = signal(SIGINT  , alya2signal_sigint  , flag )
    status = signal(SIGTERM , alya2signal_sigterm , flag ) 
    status = signal(SIGUSR1  , alya2signal_sigusr  , flag )
    status = signal(SIGUSR2  , alya2signal_sigusr  , flag )
#elif defined __ibmxl__
    call signal(SIGINT  , alya2signal_sigint  )     
    call signal(SIGTERM , alya2signal_sigterm )   
#elif defined __GNUC__
    status = signal(SIGINT  , alya2signal_sigint  )
    status = signal(SIGTERM , alya2signal_sigterm )
    status = signal(SIGUSR1 , alya2signal_sigusr )
    status = signal(SIGUSR2 , alya2signal_sigusr )
#endif
#endif
    
  end subroutine alya2signal

  subroutine alya2signal_sigint()

    print *,'SIGNAL SIGINT'
    kfl_preli = 1
    !kfl_timei = 0
    kfl_stop  = 1
    !call Restar(ITASK_WRITE_RESTART)
    !call runend('O.K.!')
    !print*,'SIGNAL SIGINT'

  end subroutine alya2signal_sigint

  subroutine alya2signal_sigterm()

    print *,'SIGNAL SIGTERM'
    kfl_preli = 1
    !kfl_timei = 0
    kfl_stop  = 1
    !call Restar(ITASK_WRITE_RESTART)
    !call runend('O.K.!')

  end subroutine alya2signal_sigterm

  subroutine alya2signal_sigusr()
    print *,'SIGNAL SIGUSR: Alya will stop after current time-step, please wait'
    kfl_preli = 1
    !kfl_timei = 0
    kfl_stop  = 1

  end subroutine alya2signal_sigusr
end module mod_alya2signal
!> @}
