subroutine openfi(itask)
  !------------------------------------------------------------------------
  !****f* master/openfi
  ! NAME 
  !    openfi
  ! DESCRIPTION
  !    This subroutine gets ALL the file names to be used by Alya in two
  !    possible ways, opens or closes them:
  ! 
  !    1. Recalling them from the environment, when Alya is launched
  !    encapsulated in a shell script, or
  ! 
  !    2. Composing the names out of the problem name which is given as
  !     argument when the binary file Alya is launched "naked". 
  !
  ! OUTPUT
  ! USES
  ! USED BY
  !    Turnon
  !    Turnof
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_postpr
  use mod_iofile
  use mod_memory,        only : lun_memor
  use mod_memory,        only : lun_varcount
  use mod_memory,        only : kfl_memor
  use mod_memory,        only : kfl_varcount
  use mod_ker_detection, only : ker_detection_open_file
  use mod_ker_timeline,  only : ker_timeline_open_file
  use mod_getopt,        only : getopt_t, getopt_long, longoption, optarg, getopt_long_help
  use mod_getopt,        only : kfl_check_data_file
  use mod_outfor,        only : outfor
  use mod_messages,      only : livinf
  implicit none

  integer(ip), intent(in) :: itask
  character(150),save     :: fil_outpu,fil_memor,fil_livei,fil_commu
  character(150),save     :: fil_latex,fil_gnupl,fil_syste
  character(150)          :: fil_memory
  character(150)          :: fil_varcount
  integer                 :: np
  character               :: option
  ! Set up the longopts struct to define the valid options:
  ! short option, long option, argument (0/1), short description:
  type(getopt_t) :: longopts(5) = [ &
       getopt_t('c', 'check',   0, 'Check data file'),    &
       getopt_t('f', 'file',    1, 'Specify input file'), &
       getopt_t('n', 'name',    1, 'Dirichlet/Neumann with PLE'), &
       getopt_t('h', 'help',    0, 'Print help'),         &
       getopt_t('',  'ignore',  0, '')                    ]

  if( itask == 1 ) then
     !
     ! Get options
     !
     ! scan all the command-line parameters
     ! getopt_long() returns a single character" ">","!",".", or the short-option character (e.g. "a" for -a).
     !   It also sets two 'global' variables through the SUFR_getopt module:
     !   - longOption:  the full option (e.g. "-a" or "--all") including the dashes
     !   - optArg:      the argument following the option (if required and present)
     if( INOTSLAVE ) then
        np = 0
        kfl_check_data_file = 0
        do
           option = getopt_long(longopts)

           ! Do different things depending on the option returned:
           select case(option)
           case('>')  ! Last parameter
              if(command_argument_count()==0) call outhel()  ! No parameters found - print help
              exit
           case('f')
              namda = trim(optarg)
              kfl_naked = 1
           case('c')
              kfl_check_data_file = 1
           case('n')
              continue
           case('h')
              !call getopt_long_help(longopts)
              call outhel()
           case('.')  ! Parameter is not an option (i.e., it doesn't start with "-" or "--")
              namda = trim(optarg)
              kfl_naked = 1
              np = np + 1
           case default
              select case(longoption)
              case('--ignore')  ! Note that --ignore was not given a short equivalent
              case default
              end select
           end select
        end do
        !
        ! Open input data file: Get file names
        !
        !call GETARG(one4,namda)
        !if(len(trim(namda))>0) then
        !   if(      trim(namda)=='--h'.or.trim(namda)=='--help'&
        !        .or.trim(namda)=='-h' .or.trim(namda)=='-help') then
        !      call outhel
        !   else
        !      kfl_naked=1
        !   end if
        !else
        !   kfl_naked=0
        !   call GET_ENVIRONMENT_VARIABLE('ALYA_NAME',namda)
        !   namda=trim(namda)
        !end if
        !
        ! Open data file
        !     
        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR011',fil_pdata)
        else
           if(adjustl(trim(namda))=='connard') then
              call runend('CONNARD TOI-MEME')
           end if
           fil_pdata = adjustl(trim(namda))//'.dat'
        end if
        call iofile(zero,lun_pdata,fil_pdata,'DATA','old')
     end if

  else if( itask == 2 ) then
     !
     ! Open output files
     !
     if( INOTSLAVE ) then
        if ( kfl_naked == 0 ) then
           call get_environment_variable('FOR012',fil_outpu) 
           call GET_ENVIRONMENT_VARIABLE('FOR013',fil_memor) 
           call GET_ENVIRONMENT_VARIABLE('FOR014',fil_conve) 
           call GET_ENVIRONMENT_VARIABLE('FOR015',fil_postp) 
           call GET_ENVIRONMENT_VARIABLE('FOR016',fil_livei)    
           call GET_ENVIRONMENT_VARIABLE('FOR017',fil_rstar)    
           call GET_ENVIRONMENT_VARIABLE('FOR018',fil_latex)    
           call GET_ENVIRONMENT_VARIABLE('FOR019',fil_gnupl)
           call GET_ENVIRONMENT_VARIABLE('FOR020',fil_commu)
           call GET_ENVIRONMENT_VARIABLE('FOR024',fil_binar)
           call GET_ENVIRONMENT_VARIABLE('FOR028',fil_syste)
           call GET_ENVIRONMENT_VARIABLE('FOR035',fil_quali) 
           call GET_ENVIRONMENT_VARIABLE('FOR045',fil_rstib)
           call GET_ENVIRONMENT_VARIABLE('FOR034',fil_memory)
           !call GET_ENVIRONMENT_VARIABLE('FOR050',fil_time)    
           !call GET_ENVIRONMENT_VARIABLE('FOR033',fil_detec)    
           !call GET_ENVIRONMENT_VARIABLE('FOR039',fil_rstla) 
           !call GET_ENVIRONMENT_VARIABLE('FOR040',fil_posla)
           !call GET_ENVIRONMENT_VARIABLE('FOR041',fil_cvgla)
        else if ( kfl_naked == 1 ) then
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
           fil_pospl  = adjustl(trim(namda))                      ! Unit?
           fil_rstib  = adjustl(trim(namda))//'-IB.rst'           ! Unit?
           fil_memory = adjustl(trim(namda))//'-memory.res'       ! Unit 34
           !fil_time  = adjustl(trim(namda))//'.tim'              ! Unit 50
           !fil_rstla = adjustl(trim(namda))//'-lagrangian.rst'
           !fil_posla = adjustl(trim(namda))//'-lagrangian.res'
           !fil_cvgla = adjustl(trim(namda))//'-lagrangian.cvg'
        end if
        !
        ! Open permanent files
        !
        if( kfl_rstar == 2 ) then
           call iofile(zero,lun_outpu,fil_outpu,'RUN EVOLUTION',      'old','formatted','append')
           call iofile(zero,lun_conve,fil_conve,'CONVERGENCE HISTORY','old','formatted','append')
           !call iofile(zero,lun_time,fil_time,'TIME HISTORY','old','formatted','append')
           call iofile(zero,lun_syste,fil_syste,'SYSTEM INFO','old','formatted','append')
           call iofile(zero,lun_memory,fil_memory,'MEMORY EVOLUTION','old','formatted','append')
        else
           call iofile(zero,lun_outpu,fil_outpu,'RUN EVOLUTION')
           call iofile(zero,lun_conve,fil_conve,'CONVERGENCE HISTORY')
           !call iofile(zero,lun_time,fil_time,'TIME HISTORY')
           call iofile(zero,lun_syste,fil_syste,'SYSTEM INFO')
           call iofile(zero,lun_memory,fil_memory,'MEMORY EVOLUTION')
        end if
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
           call livinf(-1_ip,' ',zero)
           call livinf(-2_ip,' ',zero)
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
