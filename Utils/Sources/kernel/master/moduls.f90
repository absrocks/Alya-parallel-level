subroutine moduls(jtask)
  !-----------------------------------------------------------------------
  !****f* master/moduls
  ! NAME
  !    moduls
  ! DESCRIPTION
  !    This routine calls the modules
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use def_inpout
  use mod_ker_proper
  use mod_ker_timeline
  use mod_parall,          only :  par_code_zone_subd_to_color
  use mod_parall,          only :  PAR_COMM_COLOR_ARRAY
  use mod_parall,          only :  PAR_COMM_WORLD,commd
  use mod_parall,          only :  PAR_COMM_MY_CODE_ARRAY
  use mod_parall,          only :  PAR_COMM_MY_CODE4
  use mod_communications,  only :  PAR_BARRIER
  use mod_coupling_driver, only :  COU_DRIVER
  use mod_messages,        only :  livinf
  use mod_timings,         only :  timings_doiter
  implicit none
  integer(ip), intent(in)  :: jtask
  integer(ip)              :: iorde,pmodu,itask,kmodu
  real(rp)                 :: time1,time2,time3,time4,time5
  integer(ip)              :: icolo

  ITASK_CURREN = abs(jtask)

  iorde = 0
  itask = abs(jtask)
  if( jtask >= 0 ) then
     pmodu = mmodu
  else
     pmodu = 1
  end if
  
  do while( iorde < pmodu .and. kfl_reset /= 1)
     iorde = iorde + 1

     if( itask == ITASK_TURNON ) then
        if( jtask > 0 ) modul = iorde
        if( kfl_modul(modul) == 1 ) then
           if( modul > 0 ) then
              if( kfl_itask(itask,modul) == 0 ) then
                 if( modul /= mmodu ) call livinf(51_ip,' ',modul)
              end if
           end if
           lispa = 0
           lisda = momod(modul) % lun_pdata ! Reading file
           lisre = momod(modul) % lun_outpu ! Writing file
        else
           modul = -1
        end if

     else if( itask == ITASK_DOITER .or. &
          &   itask == ITASK_CONCOU .or. &
          &   itask == ITASK_INIUNK .or. &
          &   itask == ITASK_OUTPUT .or. &
          &   itask == ITASK_ENDSTE ) then

        if( jtask > 0 ) modul = lmord(iorde,iblok)

     else if( itask == ITASK_FILTER ) then

        if( jtask > 0 ) then
           modul = modul
           iorde = mmodu
        end if

     else

        if( jtask > 0 ) then
           if( kfl_modul(iorde) == 1 ) then
              modul = iorde
           else
              modul = -1
           end if
        end if

     end if
     !
     ! Delayed module
     !
     if(  itask == ITASK_DOITER .or. &
          itask == ITASK_ENDSTE .or. &
          itask == ITASK_CONCOU ) then

        if( modul > 0 ) then
           if( kfl_delay(modul) /= 0 ) then
              modul = -1
           else
              continue
           end if
        end if

     end if
     !
     ! Pointers
     !
     if( kfl_timin == 1 ) call Parall(20_ip)
     call cputim(time1)

     if( modul > 0 ) then

        cpu_solve =  0.0_rp
        cpu_eigen =  0.0_rp
        call moddef(9_ip)
        
        if( itask == ITASK_TIMSTE ) then
           if( ittim >= ndela(modul) ) kfl_delay(modul) = 0
        end if

        if( kfl_modul(modul) == 1 ) then
           if( itask == ITASK_INIUNK ) then
              if( kfl_itask(itask,modul) == 0 ) call livinf(77_ip,' ',0_ip)
           else if( itask == ITASK_DOITER ) then
              !if( kfl_itask(itask,modul) == 0 ) call livinf(15_ip,' ',modul)
           end if
        else if( kfl_modul(modul) == 0 ) then
           modul = -1
        end if
        !
        ! Beginning of inner iterations
        !
        if( itask == ITASK_BEGITE ) then
           itinn(modul)  = 0
        end if
        !
        ! Solve module only at beginning
        !
        if( modul > 0 ) then
           if( kfl_solve(modul) == AT_BEGINNING .and. jtask > 0 ) then
              if(  itask == ITASK_TIMSTE .or. &
                   itask == ITASK_BEGSTE .or. &
                   itask == ITASK_DOITER .or. &
                   itask == ITASK_CONCOU .or. &
                   itask == ITASK_CONBLK .or. &
                   itask == ITASK_NEWMSH .or. &
                   itask == ITASK_ENDSTE ) then
                 modul = -1
              end if
           end if
        end if
        !
        ! Task can only be carried out once
        !
        if( modul > 0 ) then
           if(  itask == ITASK_REAPRO .or. &
                itask == ITASK_TURNON .or. &
                itask == ITASK_INIUNK ) then
              if( itask < 1 .or. itask > 20 ) then
                 call runend('MODULS: WRONG ITASK')
              end if
              if( kfl_itask(itask,modul) == 1 ) then
                 modul = -1
              else
                 kfl_itask(itask,modul) = 1
              end if
           end if
        end if

     end if
     !
     ! Do not recompute geometrical arrays
     !
     kfl_domar = 0
     kmodu     = modul
     !
     ! Point to zonal communication arrays and use zonal MPI communcator PAR_COMM_WORLD
     ! If module should not be called, put MODUL=-2
     !
     if( modul > 0 .and. itask == ITASK_DOITER .and. IPARALL .and. nzone > 1 ) then
        current_zone = lzone(modul)
        if( I_AM_IN_ZONE(current_zone) ) then
           icolo             =  par_code_zone_subd_to_color(current_code,current_zone,0_ip)
           commd             => PAR_COMM_COLOR_ARRAY(icolo)
           PAR_COMM_WORLD    =  PAR_COMM_COLOR_ARRAY(icolo) % PAR_COMM_WORLD
           PAR_COMM_MY_CODE4 =  int(PAR_COMM_WORLD,4)
        else
           modul = -2
        end if
     end if
     !
     ! Coupling before calling a module
     !
     if( modul >= 1 ) call COU_DRIVER(ITASK_BEFORE,itask)
     !
     ! Timing for output
     !
     if( modul >= 1 .and. itask == ITASK_OUTPUT ) then
        !call PAR_BARRIER()
        call cputim(time4)
     end if
     !
     ! Call module
     !
     lastm_ker = modul  ! module being solved for ker_proper
     !
     ! Timeline
     !
     if( kfl_timeline /= 0 .and. modul > 0 .and. modul < ID_KERMOD .and. itask >= ITASK_TURNON .and. itask <= ITASK_TURNOF ) then
        call ker_timeline(0_ip,'INI_'//trim(task_name(itask))//'_'//trim(namod(modul)),iblok)
     end if

     select case (modul)
     case (-2_ip)
        continue
     case (-1_ip)
        continue
     case ( 0_ip)
        goto 10
     case (ID_NASTIN)
        call Nastin(itask) ! NASTIN
     case (ID_TEMPER)
        call Temper(itask) ! TEMPER
     case (ID_TURBUL)
        call Turbul(itask) ! TURBUL
     case (ID_EXMEDI)
        call Exmedi(itask) ! EXMEDI
     case (ID_NASTAL)
        call Nastal(itask) ! NASTAL
     case (ID_ALEFOR)
        call Alefor(itask) ! ALEFOR
     case (ID_LATBOL)
        call Latbol(itask) ! LATBOL
     case (ID_SOLIDZ)
        call Solidz(itask) ! SOLIDZ
     case (ID_GOTITA)
        call Gotita(itask) ! GOTITA
     case (ID_WAVEQU)
        call Wavequ(itask) ! WAVEQU
     case (ID_LEVELS)
        call Levels(itask) ! LEVELS
     case (ID_QUANTY)
        call Quanty(itask) ! QUANTY
     case (ID_PARTIS)
        call Partis(itask) ! PARTIS
     case (ID_CHEMIC)
        call Chemic(itask) ! CHEMIC
     case (ID_HELMOZ)
        call Helmoz(itask) ! HELMOZ
     case (ID_IMMBOU)
        call Immbou(itask) ! IMMBOU
     case (ID_RADIAT)
        call Radiat(itask) ! RADIAT
     case (ID_CASIMI)
        call Casimi(itask) ! CASIMI
     case (ID_POROUS)
        call Porous(itask) ! POROUS
     case (ID_SOLFE2)
        call Solfe2(itask) ! SOLFE2
     !case (ID_XXXXXX)
     !   call Xxxxxx(itask) ! XXXXXX
     case (ID_NEUTRO)
        call Neutro(itask) ! NEUTRO
     case (ID_KERMOD)
        call Kermod(itask) ! KERMOD
        !if(itask==Itask_TURNON) call runend('O.K.!')
     case (ID_INSITU)
        call Insitu(itask) ! INSITU
     case default
        print *,'Unknown module called.'
     end select
     !
     ! Timeline
     !
     if( kfl_timeline /= 0 .and. modul > 0 .and. modul < ID_KERMOD .and. itask >= ITASK_TURNON .and. itask <= ITASK_TURNOF ) then
        call ker_timeline(0_ip,'END_'//trim(task_name(itask))//'_'//trim(namod(modul)),itinn(modul))
     end if
     !
     ! Coupling after calling a module
     !
     if( modul >= 1 ) call COU_DRIVER(ITASK_AFTER,itask)
     !
     ! Timing for output
     !
     if( modul >= 1 .and. itask == ITASK_OUTPUT ) then
        call cputim(time5)
        cpu_modul(ITASK_OUTPUT,modul) = cpu_modul(ITASK_OUTPUT,modul) + time5 - time4
     end if
     !
     ! Recover global communication arrays and PAR_COMM_MYCODE
     !
     if( itask == ITASK_DOITER .and. IPARALL .and. nzone > 1 ) then
        if( I_AM_IN_ZONE(current_zone) ) then
           icolo             =  par_code_zone_subd_to_color(current_code,current_zone,0_ip)
           commd             => PAR_COMM_MY_CODE_ARRAY(1)
           PAR_COMM_WORLD    =  PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD
           PAR_COMM_MY_CODE4 =  int(PAR_COMM_WORLD,4)
        end if
     end if
     !
     ! Update geometrical arrays if a module has changed them
     !
     if( ( nzone_par > 1 .and. kmodu == ID_ALEFOR ) .or. kmodu == ID_IMMBOU ) then
        call parari('MAX',0_ip,1_ip,kfl_domar)
     end if
     if( kfl_domar == 1 ) then
        call domarr(2_ip)
        kfl_domar = 0
     end if
     !
     ! Update properties
     !
     if( itask == ITASK_DOITER .or. itask == ITASK_ENDSTE ) then
        call ker_updpro(itask)
     else if (itask == ITASK_INIUNK ) then
        call ker_updpro() ! forces update (only thius time)
     end if
     !
     ! Timings
     !
     if( modul > 0 ) then
        if( kfl_modul(modul) /= 0 ) then
           call cputim(time2)
           time3                             = time2 - time1
           cpu_modul(itask,modul)            = cpu_modul(itask,modul)               + time3      ! Current task
           cpu_modul(CPU_TOTAL_MODULE,modul) = cpu_modul(CPU_TOTAL_MODULE,modul)    + time3      ! Total
           if( ITASK_CURREN == ITASK_DOITER ) then
              cpu_modul(CPU_SOLVER,modul)       = cpu_modul(CPU_SOLVER,modul)       + cpu_solve  ! Algebraic solver
              cpu_modul(CPU_EIGEN_SOLVER,modul) = cpu_modul(CPU_EIGEN_SOLVER,modul) + cpu_eigen  ! Eigen solver
              cpu_modul(CPU_DOITER,modul)       = time3
              call timings_doiter(modul)
           end if
        end if
     end if
  end do
  !
  ! Current module is kernel!
  !
10 continue
  modul =  0
  call moddef(9_ip)

end subroutine moduls
