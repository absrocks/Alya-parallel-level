!-----------------------------------------------------------------------
!> @addtogroup Master
!> @{
!> @file    openfi.f90
!> @author  houzeaux
!> @date    2020-05-11
!> @brief   Open file
!> @details This subroutine gets ALL the file names to be used by Alya in two
!>          possible ways, opens or closes them:
!>     
!>          1. Recalling them from the environment, when Alya is launched
!>          encapsulated in a shell script, or
!>     
!>          2. Composing the names out of the problem name which is given as
!>          argument when the binary file Alya is launched "naked". 
!>          
!> @} 
!-----------------------------------------------------------------------

subroutine openfi(itask)

  use def_parame
  use def_master
  use def_domain
  use def_postpr
  use mod_iofile
  use def_mpio
  use mod_memory,        only : lun_memor
  use mod_memory,        only : lun_varcount
  use mod_memory,        only : kfl_memor
  use mod_memory,        only : kfl_varcount
  use mod_ker_detection, only : ker_detection_open_file
  use mod_ker_timeline,  only : ker_timeline_open_file
  use mod_get_options,   only : get_options_help
  use mod_get_options,   only : get_options_this
  use mod_get_options,   only : getopt_t
  use mod_outfor,        only : outfor
  use mod_messages,      only : livinf
  use mod_strings,       only : string_to_integer
  use mod_messages,      only : messages_live
  implicit none

  integer(ip),    intent(in) :: itask
  character(7)               :: statu
  character(11)              :: forma
  character(6)               :: posit
  character(150), save       :: fil_outpu,fil_memor,fil_livei,fil_commu
  character(150), save       :: fil_latex,fil_gnupl,fil_syste
  character(150)             :: fil_memory
  character(150)             :: fil_varcount
  integer(ip)                :: istat
  integer(4)                 :: narg,iarg
  character(99)              :: arg
  character(99)              :: optarg
  character(1)               :: arg1
  type(getopt_t)             :: thisopt
  type(getopt_t)             :: longopts(11) = [ &
       getopt_t('k',  'check',         0_4, 'Check data file'),                    &
       getopt_t('f',  'file',          1_4, 'Specify input file'),                 &
       getopt_t('n',  'name',          1_4, 'Dirichlet/Neumann with PLE'),         &
       getopt_t('h',  'help',          0_4, 'Print help'),                         &
       getopt_t('r',  'read',          0_4, 'Read partition'),                     &
       getopt_t('w',  'write',         0_4, 'Write partition'),                    &
       getopt_t('e',  'export',        0_4, 'Export geometry in MPIO format') ,    &
       getopt_t('c',  'read-rst',      0_4, 'Read restart to continue the run'),   &
       getopt_t('i',  'read-rst-init', 0_4, 'Read restart to initialize the run'), &
       getopt_t('p',  'write-rst',     0_4, 'Write restart') ,                     &
       getopt_t('t',  'time-steps',    1_4, 'Number of time steps')                ]

  if( itask == 1 .and. INOTSLAVE ) then
     !
     ! Open data file
     !
     narg = command_argument_count()
     iarg = 0
     if( narg == 0 ) call get_options_help(longopts)

     do while( iarg < narg )
        iarg = iarg + 1
        call get_command_argument(iarg, arg)

        thisopt = get_options_this(arg,longopts)

        if( thisopt % reqarg == 1 ) then
           iarg = iarg + 1
           call get_command_argument(iarg, optarg)
        end if

        select case ( thisopt % short )

        case ('f')
           !
           ! Problem name
           !
           namda = trim(optarg)

        case ('t')
           !
           ! Number of time steps
           !
           mitim = string_to_integer(optarg,istat)
           if( istat /= 0 ) call runend('OPENFI: WRONG NUMBER OF TIME STEPS IN COMMAND LINE')

        case ('k')
           !
           ! Check data file
           !
           kfl_check_data_file = 1
           
        case ('n')

           continue

        case ('h')
           !
           ! Help
           !
           call get_options_help(longopts)

        case ('p')
           !
           ! Write restart files: preliminary run
           !
           kfl_preli = 1
           nprit     = huge(1_ip)

        case ('c')
           !
           ! Write restart files (continue run)
           !
           kfl_rstar = 2

        case ('i')
           !
           ! Write restart files (initial condition)
           !
           kfl_rstar = 1

        case ('w')
           !
           ! Write restart
           !
           kfl_ptask = 0

        case ('r')
           !
           ! Read restart
           !
           kfl_ptask = 2

        case ('e')
           !
           ! Export mesh in MPIO format
           !
           mpio_flag_geometry        = PAR_MPIO_OFF
           mpio_flag_geometry_export = PAR_MPIO_ON
           mpio_flag_post_merge      = PAR_MPIO_ON

        case default

           namda =  trim(arg)

        end select

     end do
     !
     ! Open data file
     !     
     if(adjustl(trim(namda))=='connard') then
        call runend('CONNARD TOI-MEME')
     end if
     fil_pdata = adjustl(trim(namda))//'.dat'
     if( .not. iofile_file_exists(fil_pdata) ) then
        write(6,'(a)') ''
        write(6,'(a)') 'OPENFI: WRONG PROBLEM NAME, FILE '//trim(fil_pdata)//' DOES NOT EXIST'
        write(6,'(a)') ''
        call par_finali(0_ip) 
     else
        call iofile(zero,lun_pdata,fil_pdata,'DATA','old')
     end if
     
  else if( itask == 2 ) then
     !
     ! Open output files
     !
     if( INOTSLAVE ) then
        fil_outpu  = adjustl(trim(namda))//'.log'              ! Unit 12
        fil_memor  = adjustl(trim(namda))//'.mem'              ! Unit 13
        fil_conve  = adjustl(trim(namda))//'.cvg'              ! Unit 14
        fil_livei  = adjustl(trim(namda))//'.liv'              ! Unit 16
        fil_postp  = adjustl(trim(namda))                      ! Unit 15
        fil_rstar  = adjustl(trim(namda))//'.rst'              ! Unit 17
        fil_latex  = adjustl(trim(namda))//'-latex.tex'        ! Unit 18
        fil_gnupl  = adjustl(trim(namda))//'-latex.plt'        ! Unit 19
        fil_commu  = adjustl(trim(namda))//'.com'              ! Unit 20
        fil_binar  = adjustl(trim(namda))//'.dom.bin'          ! Unit 24
        fil_syste  = adjustl(trim(namda))//'-system.log'       ! Unit 28
        fil_quali  = adjustl(trim(namda))//'-mesh-quality.res' ! Unit 35
        fil_rstib  = adjustl(trim(namda))//'-IB.rst'           ! Unit 45
        fil_memory = adjustl(trim(namda))//'-memory.res'       ! Unit 34
        !
        ! Check if restart file exists
        !
        if( kfl_rstar /= 0 .and. .not. iofile_file_exists(fil_rstar) ) then
           kfl_rstar = 0
           call messages_live('THIS IS NOT A RESTART','WARNING')
           call messages_live('THIS IS NOT A RESTART','WARNING')
           call messages_live('THIS IS NOT A RESTART','WARNING')
           call messages_live('THIS IS NOT A RESTART','WARNING')
           call messages_live('THIS IS NOT A RESTART','WARNING')
           call messages_live('KERNEL RESTART FILE DOES NOT EXIST... CONTINUING JUST IN CASE','WARNING')
        end if
        !
        ! Open permanent files
        !
        if( kfl_rstar == 2 ) then
           call iofile_restart_run(statu,forma,posit)
        else
           call iofile_normal_run(statu,forma,posit)
        end if
        call iofile(zero,lun_outpu,fil_outpu,  'RUN EVOLUTION'      ,statu,forma,posit)
        call iofile(zero,lun_conve,fil_conve,  'CONVERGENCE HISTORY',statu,forma,posit)
        call iofile(zero,lun_syste,fil_syste,  'SYSTEM INFO'        ,statu,forma,posit)
        call iofile(zero,lun_memory,fil_memory,'MEMORY EVOLUTION'   ,statu,forma,posit)
        !
        ! Open Postprocess file
        !
        if( kfl_postp_par == 1 ) then
           !
           ! Compose mesh and result file names
           !
           call outres()

        end if
        !
        ! Open live file
        !
        if(lun_livei==16) then
           if(kfl_rstar==2) then
              call iofile(zero,lun_livei,fil_livei,'LIVE INFORMATION',   'old','formatted','append')
           else
              call iofile(zero,lun_livei,fil_livei,'LIVE INFORMATION','unknown','formatted')
           end if
        end if
        !
        ! Timeline file
        !
        call ker_timeline_open_file()
        !
        ! Open latex file
        ! 
        if( kfl_latex == 1 ) then
           call iofile(zero,lun_latex,fil_latex,'LATEX')
           call iofile(zero,lun_gnupl,fil_gnupl,'GNUPLOT')
        end if
        !
        ! Open communication-with-Alya file
        !
        if( kfl_commu == 1 ) then
           call iofile(zero,lun_commu,fil_commu,'COMMUNICATION')
           write(lun_commu,'(a)') 'WRITE HERE YOUR MESSAGE TO ALYA'
        end if
        !
        ! Write header
        !
        call outfor(16_ip,lun_outpu,' ')
        call outfor(16_ip,lun_syste,' ')
     end if

  else if( itask == 3 ) then
     if( INOTSLAVE ) then
        !
        ! Close input data file
        !
        call iofile(two,lun_pdata,' ','DATA')   
        !
        ! Event detection file
        !
        call ker_detection_open_file()
     end if

  else if( itask == 7 ) then
     !
     ! Open geometry binary file
     !
     if( INOTSLAVE ) then
        call iofile(zero,lun_binar,fil_binar,'GEOMETRY BIN','unknown','unformatted')
     end if

  else if( itask == 8 ) then
     !
     ! Close geometry binary file
     !
     if( INOTSLAVE ) then
        call iofile(two,lun_binar,' ','GEOMETRY BIN') 
     end if

  else if( itask == 13 ) then
     !
     ! Open Mesh quality file
     ! 
     if( INOTSLAVE ) then
        if( kfl_quali /= 0 ) then
           call iofile(zero,lun_quali,fil_quali,'MESH QUALITY')
        end if
     end if

  else if( itask == -13 ) then
     !
     ! Open Mesh quality file
     ! 
     if( INOTSLAVE ) then
        if( kfl_quali /= 0 ) then
           call iofile(two,lun_quali,fil_quali,'MESH QUALITY')
        end if
     end if

  else if( itask == 9 ) then
     !
     ! Open memory file
     !
     kfl_memor    = abs(kfl_memor)
     kfl_varcount = abs(kfl_varcount)
     if( ISLAVE .and. kfl_outpu_par == 0 ) then
        kfl_memor    = 0
        kfl_varcount = 0
     end if
     if( kfl_memor == 1 ) then
        if ( kfl_naked == 0 ) then
           call GET_ENVIRONMENT_VARIABLE('FOR013',fil_memor)
        else if ( kfl_naked == 1 ) then
           fil_memor     = adjustl(trim(namda))//'.mem'              ! Unit 13
        end if
        if( INOTSLAVE .or. kfl_outpu_par == 1 ) then
           if( ISLAVE ) call iofile_append_tag(fil_memor,kfl_paral)
           if( kfl_rstar == 2 ) then
              call iofile(zero,lun_memor,fil_memor,'MEMORY EVOLUTION','old','formatted','append')
           else
              call iofile(zero,lun_memor,fil_memor,'MEMORY EVOLUTION')
           end if
        end if
     end if
     if( kfl_varcount == 1 ) then
        if ( kfl_naked == 0 ) then
           call GET_ENVIRONMENT_VARIABLE('FOR050',fil_varcount)
        else if ( kfl_naked == 1 ) then
           fil_varcount  = adjustl(trim(namda))//'.varcount'          ! Unit 50
        end if
        if( INOTSLAVE .or. kfl_outpu_par == 1 ) then
           if( ISLAVE ) call iofile_append_tag(fil_varcount,kfl_paral)
           if( kfl_rstar == 2 ) then
              call iofile(zero,lun_varcount,fil_varcount,'VARIABLE MEMORY COUNTER','old','formatted','append')
           else
              call iofile(zero,lun_varcount,fil_varcount,'VARIABLE MEMORY COUNTER')
           end if
        end if
     end if

  end if

end subroutine openfi
