subroutine wav_turnon
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_turnon
  ! NAME 
  !    wav_turnon
  ! DESCRIPTION
  !    This routine performs the following tasks:
  !    - Gets file names and open them.
  !    - Read data for the temperature equation.
  !    - Write some info
  !    - Check if there are errrors
  !    - Allocate memory
  ! USES
  !    wav_openfi
  !    wav_reaphy
  !    wav_reabcs
  !    wav_reanut
  !    wav_reaous
  !    wav_outinf
  !    wav_memall
  !    wav_restar
  ! USED BY
  !    Wavequ
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_wavequ
  use      mod_iofile
  implicit none
  !
  ! Initial variables
  !
  call wav_inivar(1_ip)
  !
  ! Read the physical problem
  !
  call wav_reaphy 
  !
  ! Read the numerical treatment
  !
  call wav_reanut
  !
  ! Read the output strategy
  !
  call wav_reaous
  !
  ! Service: Parall
  !
  call wav_parall(1_ip)
  !
  ! Read the boundary conditions
  !
  call wav_reabcs
  !
  ! Service: Parall
  !
  call wav_parall(2_ip)
  !
  ! Initial variables
  !
  call wav_inivar(2_ip)
  !
  ! Write info
  !
  call wav_outinf
  !
  ! Warnings and errors
  !
  call wav_outerr
  !
  ! Allocate memory
  !
  call wav_memall

end subroutine wav_turnon
