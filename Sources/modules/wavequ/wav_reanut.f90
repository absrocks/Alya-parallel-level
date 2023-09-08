subroutine wav_reanut()
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_reanut
  ! NAME 
  !    wav_reanut
  ! DESCRIPTION
  !    This routine reads the numerical treatment for TEMPER module
  ! USES
  !    ecoute
  ! USED BY
  !    wav_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_wavequ
  use def_domain
  use mod_memchk
  use mod_ecoute, only :  ecoute
  implicit none

  if(kfl_paral<=0) then
     ! 
     !  Initializations (defaults)
     !
     kfl_tiacc_wav = 2                                ! First order time integ.
     kfl_tisch_wav = 2                                ! Leap-frog
     kfl_timet_wav = 1                                ! Explicit(1)/Implicit(2)
     kfl_subgs_wav = 0                                ! Subgrid scale
     kfl_massm_wav = 0                                ! =0: Lumped mass, =1: consistent
     neule_wav     = 0                                ! Number of Euler time steps
     safet_wav     = 1.0_rp                           ! Safety factor
     sstol_wav     = 1.0e-8                           ! Steady state tolerance
     nebet_wav     = 0.25_rp                          ! Beta Newmark
     negam_wav     = 0.25_rp                          ! Gamma Newmark
     solve_sol     => solve    
     solve(1)%kfl_algso = -2                      ! Explicit scheme
     !
     ! Reach the section
     !
     call ecoute('wav_reanut')
     do while(words(1)/='NUMER')
        call ecoute('wav_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('wav_reanut')
        if(words(1)=='TIMET') then
           if(words(2)=='EXPLI') then
              kfl_timet_wav=1
           else if(words(2)=='IMPLI') then
              kfl_timet_wav=2
           end if
        else if(words(1)=='TIMEI') then
           if(exists('LEAPF')) then
              kfl_tisch_wav=2
           else if(exists('NEWMA')) then
              kfl_tisch_wav=3
              nebet_wav = getrea('BETA ',0.25_rp,'# beta parameter value for newmark method')
              negam_wav = getrea('GAMMA',0.5_rp,'# gamma parameter value for newmark method') 
              if(exists('CLOSE').or.exists('LUMPE')) then
                 kfl_massm_wav = 0
              else if(exists('CONSI')) then
                 kfl_massm_wav = 1
              end if
           end if
        else if(words(1)=='TIMEA') then
           kfl_tiacc_wav = int(param(1))
           neule_wav = getint('EULER',0_ip,'#EULER TIME STEPS')
        else if(words(1)=='SUBGR') then
           if(exists('ON   ')) kfl_subgs_wav=1
        else if(words(1)=='SAFET') then
           safet_wav = param(1)
        else if(words(1)=='STEAD') then
           sstol_wav = param(1)
        else if(words(1)=='MAXIM') then
           miinn_wav = int(param(1))
        else if(words(1)=='CONVE') then
           cotol_wav = param(1)
        else if(words(1)=='NORMO') then
           if(exists('L2   ')) then
              kfl_normc_wav = 2
           else if(exists('L1   ')) then
              kfl_normc_wav = 1
           else if(exists('L-inf')) then
              kfl_normc_wav = 0
           end if

        else if(words(1)=='ALGEB') then
           call reasol(1_ip)

        else if(words(1)=='PRECO') then 
           call reasol(2_ip)

        else if(words(1)=='SOLVE') then
           solve(1)%miter = int(param(1))

        else if(words(1)=='TOLER') then
           solve(1)%solco = param(1)

        else if(words(1)=='KRYLO') then
           solve(1)%nkryd = int(param(1))

        end if
     end do
  end if

end subroutine wav_reanut
