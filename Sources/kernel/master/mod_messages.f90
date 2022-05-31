!------------------------------------------------------------------------
!>
!> @addgroup Output
!> Toolbox for Alya output
!> @{
!> @name    ToolBox for output messages
!> @file    mod_messages.f90
!> @author  Guillaume Houzeaux
!> @date    05/02/2018
!> @brief   ToolBox for messages
!> @details ToolBox for messages to the user
!>
!------------------------------------------------------------------------

module mod_messages

  use def_kintyp,      only : ip,rp,lg,spmat,i1p,r1p
  use def_master
  use def_inpout
  use mod_ecoute,      only : ecoute
  use mod_random,      only : random_generate_number
  use mod_maths,       only : maths_day_of_week
  use mod_iofile,      only : iofile_open_unit
  use mod_iofile,      only : iofile_available_unit
  use mod_iofile,      only : iofile_flush_unit
  use mod_ansi_colors, only : ansi_colors
  use mod_ansi_colors, only : ansi_colors_name_to_code
  use mod_ansi_colors, only : ansi_colors_code_to_name
  use mod_strings,     only : integer_to_string
  implicit none

  private

  integer(ip),  parameter :: MESSAGES_SECTION = -1
  integer(ip),  parameter :: MESSAGES_WARNING = -2
  
  character(2)            :: list_colors(-5:mmodu)   
  
  character(25)           :: wspac
  character(30)           :: whea1
  character(10)           :: whea2
  character(20)           :: my_color
  character(3)            :: my_advance
  integer(ip)             :: my_module
  
  public :: messages_initialization
  public :: messages_header
  public :: messages_general
  public :: messages_live
  public :: messages_report
  public :: messages_ecoute
  public :: livinf
  public :: par_livinf

contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/01/2018
  !> @brief   Messages initialization
  !> @details Initialize header
  !>
  !----------------------------------------------------------------------

  subroutine messages_initialization()

    wspac = '                         '
    whea1 = '--|'
    whea2 = '--| ALYA '

    list_colors                   = ''
    list_colors(MESSAGES_SECTION) = ansi_colors_name_to_code('light blue')
    list_colors(MESSAGES_WARNING) = ansi_colors_name_to_code('light red')   
    !list_colors(ID_NASTIN)        = '44'
    !list_colors(ID_CHEMIC)        = '47'
    !list_colors(ID_EXMEDI)        = '45'
    !list_colors(ID_TEMPER)        = '42'
    !list_colors(ID_TURBUL)        = '43'
    !list_colors(ID_SOLIDZ)        = '41'
    !list_colors(ID_PARTIS)        = '46'
    
  end subroutine messages_initialization

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/01/2018
  !> @brief   Read colors
  !> @details Read colors from data file
  !>
  !----------------------------------------------------------------------

  subroutine messages_ecoute()

    character(200) :: name_color
    integer(ip)    :: code_color,imodu
    character(5)   :: mname

    call ecoute('RRUDAT')
    do while( words(1) /= 'ENDCO' )
       if( words(1) == 'WARNI' ) then
          if( trim(words(2)) /= '' ) then
             name_color                    = getcha_long('WARNI','red  ','#color',len(name_color,KIND=ip))             
             list_colors(MESSAGES_WARNING) = ansi_colors_name_to_code(name_color)
          else
             code_color                    = getint('WARNI',1_ip,'#color')
             list_colors(MESSAGES_WARNING) = integer_to_string(code_color)
          end if
       else if( words(1) == 'SECTI' ) then
          if( trim(words(2)) /= '' ) then
             name_color                    = getcha_long('SECTI','red  ','#color',len(name_color,KIND=ip))             
             list_colors(MESSAGES_SECTION) = ansi_colors_name_to_code(name_color)             
          else
             code_color                    = getint('SECTI',1_ip,'#color')
             list_colors(MESSAGES_SECTION) = integer_to_string(code_color)
          end if
       else
          mname = words(1)
          imodu = module_name_to_id(mname)
          if( trim(words(2)) /= '' ) then
             name_color         = getcha_long(mname,'red  ','#color',len(name_color,KIND=ip))             
             list_colors(imodu) = ansi_colors_name_to_code(name_color)
          else
             code_color         = getint(mname,1_ip,'#color')
             list_colors(imodu) = integer_to_string(code_color)
          end if
       end if
       call ecoute('RRUDAT')
    end do

  end subroutine messages_ecoute
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/01/2018
  !> @brief   Messages initialization
  !> @details Define things for output on screen, like the skin.
  !>          If skin.alyadat is present, fill in first three lines with:
  !>          1st line: One character defining the section indent character
  !>          2nd line: Chatacters for empty lines like output header
  !>          3rd line: 10 Characters for beginning each line
  !>          Example 1:
  !>          1st line:
  !>          2nd line:--|
  !>          3rd line:--| ALYA
  !>          Example 2:
  !>          1st line:
  !>          2nd line:>
  !>          3rd line:>
  !>          Example 2:
  !>          1st line:.
  !>          2nd line:   ALYA
  !>          3rd line:   ALYA
  !>
  !----------------------------------------------------------------------

  subroutine messages_header()

    integer(ip)  :: nunit,iostat
    integer(4)   :: iostat4,nunit4
    character(1) :: ww1

    if( IMASTER ) then

       nunit  = iofile_available_unit()
       nunit4 = int(nunit,4)
       call iofile_open_unit(nunit,'skin.alyadat','SKIN FOR ALYA OUTPUT',stato='old',IOSTAT=iostat)
       if( iostat == 0 ) then
          read(nunit4,1,iostat=iostat4) ww1
          read(nunit4,1,iostat=iostat4) whea1
          read(nunit4,1,iostat=iostat4) whea2
          wspac = repeat(ww1,len(wspac))
       end if

    end if

1   format(a)

  end subroutine messages_header
 
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/01/2018
  !> @brief   Color a message
  !> @details Color a message
  !
  !----------------------------------------------------------------------

  function messages_color(str, name_color, code_color) result(out)
    
    character(len=*),           intent(in) :: str
    character(len=*), optional, intent(in) :: name_color
    character(len=*), optional, intent(in) :: code_color
    character(len=200)                     :: out

    if( kfl_color == 0 .or. lun_livei /= 6 ) then
       out = trim(str)
    else
       out = ansi_colors(str, name_color, code_color)
    end if
    
  end function messages_color
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/01/2018
  !> @brief   Some reporting with respect to the run
  !> @details Output some reporting with respect to the data of the run
  !
  !----------------------------------------------------------------------

  subroutine messages_report()

    use def_domain, only : nelem_total

    if( IMASTER .and. nelem_total / max(npart,1_ip) < 10000 .and. nelem_total > 1000000 ) &
         call messages_live('ARE YOU SURE YOU NEED SO MANY CPUs? YOUR NUMBER OF ELEMENTS PER SUBDOMAIN IS VERY LOW!','REPORT')

  end subroutine messages_report

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/01/2018
  !> @brief   Live message
  !> @details Bridge to livinf
  !
  !----------------------------------------------------------------------

  subroutine messages_output(message)

    character(*),      intent(in) :: message

    if( isect > 3*max_secti ) return
    
    if( trim(my_advance) == 'no' .or. trim(adjustl(message)) == '(' .or. trim(adjustl(message)) == ')' ) then
       !
       ! Message is the continuation of a line
       !
       if( trim(my_color) == '' ) then
          write(lun_livei,FMT='(a)',ADVANCE=trim(my_advance))  trim(message)
       else
          write(lun_livei,FMT='(a)',ADVANCE=trim(my_advance))  trim(messages_color(trim(message),NAME_COLOR=trim(my_color)))
       end if
       
    else
       !
       ! This is a new line, put ALYA header
       !
       if( trim(my_color) == '' ) then
          write(lun_livei,FMT='(a)',ADVANCE='yes') whea2//repeat(' ',max(0_ip,isect))//trim(message)
       else
          write(lun_livei,FMT='(a)',ADVANCE='yes') whea2//repeat(' ',max(0_ip,isect))//trim(messages_color(trim(message),NAME_COLOR=trim(my_color)))
       end if
       
    end if
    
  end subroutine messages_output
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/01/2018
  !> @brief   Live message
  !> @details Bridge to livinf
  !
  !----------------------------------------------------------------------

  subroutine messages_live(message,message_type,RANK4,INT_NUMBER,REAL_NUMBER,INT_ARRAY,CHAR_ARRAY,REAL_ARRAY,FMT,COLOR,ADVANCE)

    character(*),     intent(in)           :: message
    character(*),     intent(in), optional :: message_type
    integer(4),                   optional :: RANK4
    integer(ip),      intent(in), optional :: INT_NUMBER
    real(rp),         intent(in), optional :: REAL_NUMBER
    integer(ip),      intent(in), optional :: INT_ARRAY(:)
    real(rp),         intent(in), optional :: REAL_ARRAY(:)
    character(len=*), intent(in), optional :: CHAR_ARRAY(:)
    character(len=*), intent(in), optional :: FMT
    character(len=*), intent(in), optional :: COLOR
    character(len=*), intent(in), optional :: ADVANCE
    integer(ip)                            :: kfl_paral_sav
    
    if( present(RANK4) ) then
       kfl_paral_sav = kfl_paral
       if( RANK4 == 0_4 .and. lun_livei == 6 ) kfl_paral = 0
    end if

    if( present(message_type) ) then
       
       if(      trim(message_type) == 'WARNING'       ) then
          call livinf(-17_ip,trim(message),0_ip,COLOR=COLOR,ADVANCE=ADVANCE)
       else if( trim(message_type) == 'REPORT'        ) then
          call livinf(-18_ip,trim(message),0_ip,COLOR=COLOR,ADVANCE=ADVANCE)
       else if( trim(message_type) == 'START SECTION' ) then
          call livinf( -4_ip,trim(message),0_ip,COLOR=COLOR,ADVANCE=ADVANCE)
       else if( trim(message_type) == 'END SECTION'   ) then
          call livinf( -5_ip,trim(message),0_ip,COLOR=COLOR,ADVANCE=ADVANCE)
       else if( trim(message_type) == 'MODULE'        ) then
          call livinf( 59_ip,trim(message),0_ip,COLOR=COLOR,ADVANCE=ADVANCE)
       else
          call livinf(0_ip,trim(message),0_ip,COLOR=COLOR,ADVANCE=ADVANCE)
       end if
    else
       
       if( trim(message) == '---' ) then
          
          call livinf(0_ip,repeat("-",47) ,0_ip,COLOR=COLOR)
          
       else if( present(REAL_NUMBER) ) then
          
          if( isect > 3*max_secti ) return
          if( present(FMT) ) then
             write(lun_livei,FMT=FMT) whea2//repeat(' ',max(0_ip,isect))//message,REAL_NUMBER
          else
             write(lun_livei,FMT='(a,1x,e13.6)') whea2//repeat(' ',max(0_ip,isect))//message,REAL_NUMBER
          end if

          
       else if( present(INT_NUMBER) ) then
          
          if( isect > 3*max_secti ) return
          if( present(FMT) ) then
             write(lun_livei,FMT=FMT) whea2//repeat(' ',max(0_ip,isect))//message,INT_NUMBER
          else
             write(lun_livei,FMT='(a,1x,i6)') whea2//repeat(' ',max(0_ip,isect))//message,INT_NUMBER
          end if
          
       else if( present(REAL_ARRAY) ) then
          
          if( isect > 3*max_secti ) return
          if( present(FMT) ) then
             write(lun_livei,FMT=FMT) whea2//repeat(' ',max(0_ip,isect))//message,REAL_ARRAY
          else
             write(lun_livei,FMT='(a,300(1x,e13.6))') whea2//repeat(' ',max(0_ip,isect))//message,REAL_ARRAY
          end if
          
       else if( present(INT_ARRAY) ) then
          
          if( isect > 3*max_secti ) return
          if( present(FMT) ) then
             write(lun_livei,FMT=FMT) whea2//repeat(' ',max(0_ip,isect))//message,INT_ARRAY
          else
             write(lun_livei,FMT='(a,300(1x,i6))') whea2//repeat(' ',max(0_ip,isect))//message,INT_ARRAY
          end if
          
       else if( present(CHAR_ARRAY) ) then
          
          if( isect > 3*max_secti ) return
          if( present(FMT) ) then
             write(lun_livei,FMT=FMT) whea2//repeat(' ',max(0_ip,isect))//message,CHAR_ARRAY
          else
             write(lun_livei,FMT='(a,300(1x,a))') whea2//repeat(' ',max(0_ip,isect))//message,CHAR_ARRAY
          end if
          
       else
          call livinf(0_ip,trim(message),0_ip,COLOR=COLOR,ADVANCE=ADVANCE)
       end if
    end if
    
    if( present(RANK4) ) then
       kfl_paral = kfl_paral_sav
    end if
    
  end subroutine messages_live

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/01/2018
  !> @brief   Final output message
  !> @details Output message when finalizing Alya
  !
  !----------------------------------------------------------------------

  subroutine messages_general(message)

    character(*), intent(out) :: message

    character(20)             :: wuser
    character(10)             :: wtime
    integer                   :: values(8)
    integer(ip)               :: imess
    character(9)              :: wday
    logical(lg)               :: if_late,if_week_end
    character(200)            :: messa(50)
    real(rp)                  :: rr,time1

    if_late     = .false.
    if_week_end = .false.
    message     = ''
    messa       = ''
    imess       = 0

    call get_environment_variable("USER",wuser)
    if( trim(wuser) == '' ) wuser = 'User'
    !
    ! Working hours
    !
    call date_and_time(VALUES=values,TIME=wtime)
    wday = maths_day_of_week(int(values(3),8),int(values(2),8),int(values(1),8))

    if(  wtime(1:2) == '22' .or. &
         wtime(1:2) == '23' .or. &
         wtime(1:2) == '24' .or. &
         wtime(1:2) == '01' .or. &
         wtime(1:2) == '02' .or. &
         wtime(1:2) == '03' .or. &
         wtime(1:2) == '04' .or. &
         wtime(1:2) == '05' .or. &
         wtime(1:2) == '06' ) then
       if_late = .true.
    end if
    if( trim(wday) == 'Saturday' .or. trim(wday) == 'Sunday' ) then
       if_week_end = .true.
    end if
    if( if_late .and. if_week_end ) then
       imess        = imess + 1
       messa(imess) = 'You are a hard worker '//trim(wuser)//', working late in the week-end!!!'
    else if( if_week_end ) then
       imess        = imess + 1
       messa(imess) = trim(wuser)//', are you sure you should work on the week-end...?'
    else if( if_late ) then
       imess        = imess + 1
       messa(imess) = 'Hey '//trim(wuser)//', do you have a deadline to work so late?'
    end if
    !
    ! User
    !
    if( trim(wuser) == 'bsc21903' ) then
       imess        = imess + 1
       messa(imess) = 'Hola '//trim(wuser)//', debuggeando la TestSuite?'
    else if( trim(wuser) == 'houzeaux' ) then
       imess        = imess + 1
       messa(imess) = 'On bosse dur '//trim(wuser)//'?'
    else if( trim(wuser) == 'houzeaux' ) then
       imess        = imess + 1
       messa(imess) = 'On bosse dur '//trim(wuser)//'?'
    end if
    !
    ! Example
    !
    if( trim(namda(1:6)) == 'cavity' .or. trim(namda(1:6)) == 'CAVITY' .or. trim(namda(1:3)) == 'cav' ) then
       imess        = imess + 1
       messa(imess) = 'We have run so many times this cavity flow '//trim(wuser)
    end if
    if( index(namda,'spray') /= 0 .or. index(namda,'SPRAY') /= 0 ) then
       if( trim(wuser) == 'bsc21304' .or. trim(wuser) == 'bsc21808' ) then
          imess        = imess + 1
          messa(imess) = 'Dear '//trim(wuser)//', are you willing to design a spray? Maybe you should look in google... try https://www.google.com/search?q=spray'
       end if
    end if
    !
    ! Modules
    !
    if( kfl_modul(ID_PARTIS) == 1 ) then
       imess        = imess + 1
       messa(imess) = 'Hope you have not lost too many particles '//trim(wuser)//'...'//&
            & ' If this is the case, you may have strange elements!'
    end if
    if( sum(kfl_modul(1:mmodu-1)) >= 3 ) then
       imess        = imess + 1
       messa(imess) = 'Whoo! You are coupling lots of modules! Greetings if it worked ;o)'
    end if
    !
    ! Others
    !
#ifdef ALYA_DLB
    imess        = imess + 1
    messa(imess) = 'Estas usando DLB '//trim(wuser)//', quieres ahorrar tiempo?'
#endif
#ifdef ALYA_OMPSS
    imess        = imess + 1
    messa(imess) = 'Usando OmpSS '//trim(wuser)//', que atrevido!'
#endif
    if( npart > 10000 ) then
       imess        = imess + 1
       messa(imess) = trim(wuser)//', are you sure you need so many CPUs ;o) ?'
    end if
    call cputim(time1)
    time1 = time1 - cpu_initi
    if( time1 > 3600.0_rp*12.0_rp ) then
       imess        = imess + 1
       messa(imess) = 'Your run is quite long '//trim(wuser)//'!'
    end if
    !
    ! Choose a message among the generated
    !
    if( imess > 0 ) then
       rr    = random_generate_number()
       imess = (int(rr*real(imess,rp),ip) + 1_ip)
       imess = max(imess,1_ip)
       imess = min(imess,size(messa,KIND=ip))
       rr    = random_generate_number()
       if( rr > 0.8_rp ) message = messa(imess)
    end if

  end subroutine messages_general

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-09
  !> @brief   Parallelization messages
  !> @details ?Parallelization messages
  !> 
  !-----------------------------------------------------------------------
  
  subroutine par_livinf(itask,message,inume)

    use def_domain, only : nedge
    use def_domain, only : medge
    use def_parall
    use mod_parall
    integer(ip),             intent(in) :: itask
    character(*),            intent(in) :: message
    integer(ip),   optional, intent(in) :: inume
    character(300)                      :: messa,messb
    character(20)                       :: mess1,mess2

    if( INOTSLAVE .and. lun_livei /= 0 ) then     

       if(itask==0) then
          if(inume==nproc_par-1) then
             messa='ALL SLAVES HAVE INITIATED'
          else
             mess1=intost(nproc_par-inume-1_ip)
             call runend(trim(mess1)//'PROCESSES HAVE NOT BEEN INITIATED')
          end if

       else if(itask==1) then
          messa='START MESH PARTITION'
          inews = 1

       else if(itask==2) then
          messa='MASTER COMPUTES ELEMENT GRAPH'

       else if(itask==3) then
          if( kfl_partition_par == PAR_METIS4 ) then
             if(kfl_parti_par==1) then
                messa='MASTER PARTITIONS ELEMENT GRAPH WITH METIS USING NODES CONNECTIVITY'
             else
                messa='MASTER PARTITIONS ELEMENT GRAPH WITH METIS USING FACE CONNECTIVITY'
             end if
          else if( kfl_partition_par == PAR_SFC ) then
             messa='MASTER PARTITIONS MESH USING A HILBERT SPACE FILLING CURVE'
          else if( kfl_partition_par == PAR_ZOLTAN ) then
             messa='MASTER PARTITIONS MESH USING ZOLTAN A HSFC'
          else if( kfl_partition_par == PAR_ORIENTED_BIN ) then
             messa='MASTER PARTITIONS MESH USING AN ORIENTED BIN'
          else if( kfl_partition_par < 0 ) then
             messa='MASTER PARTITIONS MESH USING A FIELD'
          end if

       else if(itask==4) then
          messa='MASTER COMPUTES COMMUNICATION STRATEGY'

       else if(itask==5) then
          messa='MASTER COMPUTES PERMUTATION ARRAYS'

       else if(itask==6) then
          if(kfl_ptask==0) then
             messa='MASTER WRITES MESH AND PARTITION DATA IN RESTART FILES'
          else if(kfl_ptask==1) then
             messa='MASTER/SLAVES EXCHANGE MESH AND PARTITION DATA'
          else if(kfl_ptask==2) then
             messa='MASTER/SLAVES READ MESH AND PARTITION DATA FROM RESTART FILE'
          end if
          if(kfl_fileh_par==1.and.kfl_ptask/=1) messa=trim(messa)//' WITH HIERARCHY'

       else if(itask==7) then
          continue

       else if(itask==8) then
          messa='MASTER/SLAVES READ MESH AND PARTITION DATA FROM RESTART FILE'

       else if(itask==9) then
          messa='MASTER OUTPUT MESH'

       else if(itask==10) then
          if(inume/=0) then
             call runend('SOME SLAVES COULD NOT READ THEIR RESTART FILES')
          else
             messa='MASTER/SLAVES: ALL SLAVES HAVE READ THEIR RESTART FILES'
          end if

       else if(itask==11) then
          mess1=intost(nproc_par-1_ip)
          messa='CHECK MPI. 1 MASTER + '//trim(mess1)//' SUBDOMAINS'

       else if(itask==12) then
          messa='MPI IS WORKING WELL'

       else if(itask==13) then
          messa='MPI IS WORKING WELL'

       else if(itask==14) then
          messa=trim(namod(modul))//': MASTER SENDS TO SLAVES GROUPS OF DEFLATED CG'

       else if(itask==18) then
          messa='MASTER ORDERS INTERIOR AND BOUNDARY NODES'

       else if (itask == 20) then
          messa='PARALL PREPROCESS -->  '//trim(message)

       else if (itask == 21) then
          messa='SLAVE START PARTITION PREPROCESS -->  '//trim(message)

       else if (itask == 22) then
          messa='SLAVE ENDS PARTITION PREPROCESS -->  '//trim(message)

       else if(itask==1000) then
          messa='END MESH PARTITION'
          isect = isect - 3
       end if

       if(itask==2) then
          messb = trim(messa)
          call livinf(-2_ip,trim(messb),0_ip)

       else if(itask==1) then
          mess1 = intost(npart_par)
          messb = trim(messa)//' (# OF SUBDOMAINS= '//trim(mess1)//')'
          call livinf(0_ip,trim(messb),0_ip)

       else if(itask==7) then
          mess1 = intost(nedge)
          mess2 = intost(medge)
          messb = ' (# EDGES= '//trim(mess1)//', MAX # EDGES/ELEMENT= '//trim(mess2)//')'
          call livinf(-3_ip,trim(messb),0_ip)

       else
          
          call messages_live(trim(messa))

       end if

    end if

    !if( inews == 1 ) isect = isect + 3

1   format(' (# EDGES= ',a,', MAX # EDGES/ELEMENT= ',a,')')

  end subroutine par_livinf
  
  !-----------------------------------------------------------------------
  !> @author  Mariano Vazquez
  !> @date    16/11/1966
  !> @brief   Echo a given message on screen
  !> @details Echo a given message on screen
  !-----------------------------------------------------------------------
  
  subroutine livinf(itask,message,INT_NUMBER,REAL_NUMBER,COLOR,ADVANCE)

    use def_master
    use mod_parall, only : mapps
    use def_domain, only : nelem,nboun,npoin

    implicit none

    integer(ip),      intent(in)           :: itask
    character(len=*), intent(in)           :: message
    integer(ip),      intent(in), optional :: INT_NUMBER
    real(rp),         intent(in), optional :: REAL_NUMBER
    character(len=*), intent(in), optional :: COLOR
    character(len=*), intent(in), optional :: ADVANCE
    character(300)                         :: messa,dumml,mess1,mess2,mess3
    integer(ip)                            :: ii,inume
    real(rp)                               :: rnume
    character(1)                           :: wbyte(3)
    real(rp)                               :: rbyte
    character(6)                           :: my_fmt
    
    !-------------------------------------------------------------------
    !
    ! Optional arguments
    !
    !-------------------------------------------------------------------
    
    my_color   = optional_argument(''      , COLOR       ) ! Name (red, green...)
    my_advance = optional_argument('yes'   , ADVANCE     ) ! Adavnce (no, yes)
    inume      = optional_argument(0_ip    , INT_NUMBER  )
    rnume      = optional_argument(0.0_rp  , REAL_NUMBER )

    !-------------------------------------------------------------------
    !
    ! Impose module color if given by user
    !
    !-------------------------------------------------------------------
    
   if( modul > 0 .and. modul < mmodu .and. .not. present(COLOR) ) then
      my_color = ansi_colors_code_to_name(list_colors(modul))
    end if
    
    !-------------------------------------------------------------------
    !
    ! Form message
    !
    !-------------------------------------------------------------------

    if( kfl_paral <= 0 .and. lun_livei /= 0 ) then

       if( itask == 1 ) then
          !
          ! Starting Alya, first message!
          !
          if(kfl_rstar == 0 ) then
             write(lun_livei,FMT=3,ADVANCE=my_advance) whea1
             write(lun_livei,FMT=3,ADVANCE=my_advance) whea2//'START ALYA FOR PROBLEM: '//trim(title)
             write(lun_livei,FMT=3,ADVANCE=my_advance) whea1
          else
             write(lun_livei,FMT=3,ADVANCE=my_advance) whea1
             write(lun_livei,FMT=3,ADVANCE=my_advance) whea2//'RESTART ALYA FOR PROBLEM: '//trim(title)
             write(lun_livei,FMT=3,ADVANCE=my_advance) whea1
          end if
          if( mapps > 1 ) then
             write(lun_livei,FMT=3,ADVANCE=my_advance) 'RUNNING ALYA WITH ANOTHER CODE'
             write(lun_livei,FMT=3,ADVANCE='no') 'NAMES OF CODES: '
             do ii = 1,mapps
                if( ii == mapps ) then
                   write(lun_livei,FMT=3,ADVANCE=my_advance) trim(application_names(ii))
                else
                   write(lun_livei,FMT=3,ADVANCE=my_advance) trim(application_names(ii))//', '
                end if
             end do
          else
             continue
          end if
          call iofile_flush_unit(lun_livei)
          return
          
       else if( itask == 4 ) then
          !
          ! START TIME STEP
          !
          if( kfl_timco /= -1 ) then
             messa = trim(messages_color('START TIME STEP '//trim(intost(ittim+1_ip)),CODE_COLOR=list_colors(MESSAGES_SECTION)))
             inews = 1
          else
             return
          end if
          
       else if( itask == 5 ) then          
          !
          ! START BLOCK
          !
          if(nblok>1.and.itcou==1 ) then
             messa = 'START BLOCK '//trim(intost(iblok))
             inews = 1
          else
             return
          end if
          
       else if( itask == 6 ) then
          !
          ! START COUPLING
          !
          if(micou(iblok)>1 ) then
             messa = 'START COUPLING ITERATION '//trim(intost(itcou))
             inews = 1
          else
             return
          end if
          
       else if( itask == 7 ) then
          !
          ! END COUPLING
          !
          if(micou(iblok)>1 ) then
             messa = 'END COUPLING ITERATION '
             isect = isect - 3
          else
             return
          end if
          
       else if( itask == 8 ) then
          !
          ! END BLOCK
          !
         if(nblok>1 ) then
             messa = 'END BLOCK'
             isect = isect - 3
          else
             return
          end if
          
       else if( itask == 9 ) then
          !
          ! END TIME STEP
          !
          if(kfl_timco/=-1 ) then
             messa = trim(messages_color('END TIME STEP',CODE_COLOR=list_colors(MESSAGES_SECTION)))
             isect = isect - 3
          else
             return
          end if
                                        
       else if( itask == 59 ) then
          !
          ! Add module name to the message
          !
          messa = trim(messages_color(trim(namod(modul))//': '//trim(message),CODE_COLOR=list_colors(inume)))
                    
       else if( itask == 15 .or. itask == 56 ) then
          !
          ! SOLVE MODUL
          !
          messa = trim(messages_color('SOLVE '//trim(namod(inume)),CODE_COLOR=list_colors(inume)))
          
       else if( itask == 160 ) then
          !
          ! (
          !
          messa = trim(messages_color(' (',CODE_COLOR=list_colors(modul)))

       else if( itask == 16 ) then
          !
          ! (X)
          !
          dumml = intost(inume)
          messa = trim(messages_color(' ('//trim(dumml)//')',CODE_COLOR=list_colors(modul)))
          
       else if( itask == 164 ) then
          !
          ! )
          !
          messa = trim(messages_color(')',CODE_COLOR=list_colors(modul)))
          
       else if( itask == 165 ) then
          !
          ! End of solving
          !
          messa = trim(messages_color(trim(message),CODE_COLOR=list_colors(modul)))
          
       else if( itask == 201 ) then
          
          isect = isect-3
          write(lun_livei,FMT='(a)')        whea2//repeat(' ',max(0_ip,isect))//'CRITICAL TIME OVERSHOOT'
          write(lun_livei,FMT='(a,es13.6)') whea2//repeat(' ',max(0_ip,isect))//'GOING BACK TO TIME t= ',cutim
          messa = trim('RESTARTING TIME STEP '//trim(intost(ittim)))
          isect = isect+3
          
       else if( itask == 10000 ) then
          
          write(lun_livei,FMT=3,ADVANCE=my_advance) whea1
          messa=whea1//' WARNING: '// adjustl(trim(message))
          write(lun_livei,FMT=3,ADVANCE=my_advance) whea2//repeat(' ',max(0_ip,isect))//adjustl(trim(messa))
          write(lun_livei,FMT=3,ADVANCE=my_advance) whea1
          call iofile_flush_unit(lun_livei)

       end if

       !-------------------------------------------------------------------
       !
       ! Generic messages: to be used preferently
       !
       !-------------------------------------------------------------------

       if( itask == 0  .or. itask == 42 ) then
          
          call messages_output(message)
          
       else if( itask == 97 ) then
          !
          ! Mesh multiplication
          !
          do ii = 1,3
             if( routp(ii) >= 1.0e9_rp ) then
                rbyte     = 1.0e9_rp
                wbyte(ii) = 'B'
             else if( routp(ii) >= 1.0e6_rp ) then
                rbyte     = 1.0e6_rp
                wbyte(ii) = 'M'
             else if( routp(ii) >= 1.0e3_rp ) then
                rbyte     = 1.0e3_rp
                wbyte(ii) = 'k'
             else
                rbyte     = 1.0_rp
                wbyte(ii) = ''
             end if
             routp(ii) = routp(ii) / rbyte
             ioutp(ii) = int(routp(ii),ip)
          end do

          messa = ' (NELEM= '//trim(mess1)//', NPOIN= '//trim(mess2)//', NBOUN= '//trim(mess3)//')'
          write(lun_livei,FMT=5) whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '&
               //' (NELEM= ',(routp(1)),wbyte(1) &
               //', NPOIN= ',(routp(2)),wbyte(2) &
               //', NBOUN= ',(routp(3)),wbyte(3) &
               //')'

       else if( itask == -3  ) then
          !
          ! Solidz
          !
          write(lun_livei,FMT='(a)') trim(messages_color(trim(message),CODE_COLOR=list_colors(modul))) 
          
       else if( itask == -6 ) then
          !
          ! Coupling
          !
          inews = 1
          messa = intost(inume)
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '//trim(messa)

       else if( itask == -7 ) then
          !
          ! Partis
          !
          write(lun_livei,FMT='(a,i7)') whea2//repeat(' ',max(0_ip,isect))//trim(message)//': ',inume

       else if( itask == -9 ) then
          !
          ! Partis and Chemic
          !
          if( inume < 0 ) then
             my_fmt = '(a,i7)'
          else if( inume >= 10**6 ) then
             my_fmt = '(a,i7)'
          else if( inume >= 10**5 ) then
             my_fmt = '(a,i6)' 
          else if( inume >= 10**4 ) then
             my_fmt = '(a,i5)' 
          else if( inume >= 10**3 ) then
             my_fmt = '(a,i4)' 
          else if( inume >= 10**2 ) then
             my_fmt = '(a,i3)'
          else if( inume >= 10 ) then
             my_fmt = '(a,i2)'
          else
             my_fmt = '(a,i1)'
          end if
          write(lun_livei,FMT=my_fmt) whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message),inume

       else if( itask == -13 ) then
          !
          ! Coupling
          !
          isect = isect - 3
          messa = intost(inume)
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '//trim(messa)

       else if( itask == -14 ) then
          !
          ! Coupling
          !
          messa = intost(inume)
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '//trim(messa)

       else if( itask == -16 ) then
          !
          ! Nastin
          !
          write(lun_livei,FMT='(a,e12.6)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message)//'= ',routp(1)

       else if( itask == -17  ) then
          !
          ! Warning
          !
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(messages_color('!!! WARNING: '//trim(message),CODE_COLOR=list_colors(MESSAGES_WARNING)))

       else if( itask == -18  ) then
          !
          ! Report
          !
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//'!!! REPORT: '//trim(message)

       else if( itask == -4 ) then
          !
          ! Start section
          !
          inews = 1
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(messages_color('START '//trim(message),CODE_COLOR=list_colors(MESSAGES_SECTION)))

       else if( itask == -5 ) then
          !
          ! End section
          !
          isect = isect - 3
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(messages_color('END '//trim(message),CODE_COLOR=list_colors(MESSAGES_SECTION)))

          !-------------------------------------------------------------------
          !
          ! Specific messages
          !
          !-------------------------------------------------------------------
          
       else if( itask == 4 .or. itask == 15 .or. itask == 44 .or. itask == 43 .or. itask == 56 ) then
          !
          ! Do not advance
          !
          write(lun_livei,FMT='(a)',advance='no') whea2//repeat(' ',max(0_ip,isect))//trim(messa)

       else if( itask == 165 .or. itask == 161 .or. itask == 160 .or. itask == 163 ) then
          !
          ! Number of iterations, without title
          !
          write(lun_livei,FMT='(a)',advance='no') trim(messa)

       else if( itask == 164 .or. itask == 16 ) then
          !
          ! End of iteration without title
          !
          write(lun_livei,FMT='(a)') trim(messa)

       else if( itask == 18 ) then
          !
          ! Time step
          !
          write(lun_livei,FMT='(a,es13.6)') ', t= ',cutim

       else if( itask == 98  ) then
          !
          ! Alefor imbbou
          !
          write(lun_livei,FMT='(a,e16.8E3)') whea2//repeat(' ',max(0_ip,isect))//trim(message),routp(1)

       else if( itask == 10000 ) then
          !
          ! Warning
          !
          write(lun_livei,FMT=3,ADVANCE=my_advance) whea1
          messa=whea1//' WARNING: '// adjustl(trim(message))
          write(lun_livei,FMT=3,ADVANCE=my_advance) whea2//repeat(' ',max(0_ip,isect))//adjustl(trim(messa))
          write(lun_livei,FMT=3,ADVANCE=my_advance) whea1
          call iofile_flush_unit(lun_livei)

       else if( itask == -1 ) then
          !
          ! Runend
          !
          messa = ''//trim(message)
          write(lun_livei,FMT=3,ADVANCE=my_advance) trim(messa)

       else

          write(lun_livei,FMT=3,ADVANCE=my_advance) whea2//repeat(' ',max(0_ip,isect))//trim(messa)

       end if
       
       if(itask>=0) call iofile_flush_unit(lun_livei)

    else if(kfl_paral>=-1 ) then

       if( itask == 999 ) then
          write(6,3) ' '
          messa=whea2//'ABORTED. '//trim(message)
          write(6,3) trim(messa)
       end if

    end if

    if( inews == 1 ) isect = isect + 3
    inews = 0

    return
1   format(i3,'%...')
3   format(a)
4   format(a,a11,1x,a,1x,a11,1x,a,1x,a3)
5   format(a,f5.1,a,f5.1,a,f5.1,a)

  end subroutine livinf

end module mod_messages
!> @}
