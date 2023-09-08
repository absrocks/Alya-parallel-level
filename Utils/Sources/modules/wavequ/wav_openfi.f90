subroutine wav_openfi(itask)

  !-----------------------------------------------------------------------
  !    
  ! This subroutine gets ALL the file names and open them to be used by 
  ! the module in two possible ways:
  ! 
  ! 1. Recalling them from the environment, when Alya is launched
  ! encapsulated in a shell script, or
  ! 
  ! 2. Composing the names out of the problem name which is given as argument
  ! when the binary file Alya is launched "naked". 
  !
  !-----------------------------------------------------------------------
  use      def_wavequ
  use      def_master
  use      mod_iofile
  implicit none
  integer(ip), intent(in) :: itask 
  character(150)          :: fil_pdata_wav,fil_outpu_wav
  character(150)          :: fil_conve_wav,fil_solve
  character(150)          :: fil_setse_wav,fil_setsb_wav,fil_setsn_wav
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

end subroutine wav_openfi

