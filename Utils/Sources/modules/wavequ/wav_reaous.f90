subroutine wav_reaous

  !-----------------------------------------------------------------------
  !
  ! This routine reads the output strategy for the temperature
  ! equation.
  !
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_wavequ
  use      def_domain
  use      mod_memchk
  use mod_ecoute, only :  ecoute
  implicit none

  call reaous()

end subroutine wav_reaous
    
