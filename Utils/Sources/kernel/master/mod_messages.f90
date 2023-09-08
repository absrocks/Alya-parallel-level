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

  use def_kintyp, only : ip,rp,lg,spmat,i1p,r1p
  use def_master
  use mod_random, only : random_grnd
  use mod_maths,  only : maths_day_of_week
  use mod_iofile, only : iofile_open_unit
  use mod_iofile, only : iofile_available_unit
  use mod_iofile, only : iofile_flush_unit
  implicit none

  private

  character(25) :: wspac
  character(30) :: whea1
  character(10) :: whea2

  public :: messages_initialization
  public :: messages_header
  public :: messages_general
  public :: messages_live
  public :: messages_report
  public :: livinf

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

  end subroutine messages_initialization

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

  subroutine messages_live(message,message_type,RANK4)

    character(*), intent(in)           :: message
    character(*), intent(in), optional :: message_type
    integer(4),               optional :: RANK4
    integer(ip)                        :: kfl_paral_sav

    kfl_paral_sav = kfl_paral
    if( present(RANK4) ) then
       if( RANK4 == 0_4 .and. lun_livei == 6 ) kfl_paral = 0
    end if

    if( present(message_type) ) then
       if(      trim(message_type) == 'WARNING'       ) then
          call livinf(-17_ip,trim(message),0_ip)
       else if( trim(message_type) == 'REPORT'        ) then
          call livinf(-18_ip,trim(message),0_ip)
       else if( trim(message_type) == 'START SECTION' ) then
          call livinf(-19_ip,trim(message),0_ip)
       else if( trim(message_type) == 'END SECTION'   ) then
          call livinf(-20_ip,trim(message),0_ip)
       else
          call livinf(0_ip,trim(message),0_ip)
       end if
    else
       call livinf(0_ip,trim(message),0_ip)
    end if
    kfl_paral = kfl_paral_sav

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
#ifdef __ibmxl__
    message = 'What a miracle! You succeeded in running Alya by compiling with XLF!'
#else
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
       rr    = random_grnd()
       imess = (int(rr*real(imess,rp),ip) + 1_ip)
       imess = max(imess,1_ip)
       imess = min(imess,size(messa,KIND=ip))
       rr    = random_grnd()
       if( rr > 0.8_rp ) message = messa(imess)
    end if
#endif

  end subroutine messages_general

  !-----------------------------------------------------------------------
  !> @author  Mariano Vazquez
  !> @date    16/11/1966
  !> @brief   Echo a given message on screen
  !> @details Echo a given message on screen
  !-----------------------------------------------------------------------
  subroutine livinf(itask,message,inume)

    use def_master
    use mod_parall, only : mapps
    use def_domain, only : nelem,nboun,npoin

    implicit none

    integer(ip),      intent(in) :: itask
    integer(ip),      intent(in) :: inume
    character(len=*), intent(in) :: message
    character(300)               :: messa,dumml,mess1,mess2,mess3
    integer(ip)                  :: ii
    character(1)                 :: wbyte(3)
    real(rp)                     :: rbyte

    !-------------------------------------------------------------------
    !
    ! Form message
    !
    !-------------------------------------------------------------------

    if( kfl_paral <= 0 .and. lun_livei /= 0 ) then

       if((itask==1.and.lun_livei==6).or.(itask==-1.and.lun_livei/=6) ) then
          if(kfl_rstar==0 ) then
             write(lun_livei,FMT=3) whea1
             write(lun_livei,FMT=3) whea2//'START ALYA FOR PROBLEM: '//trim(title)
             write(lun_livei,FMT=3) whea1
          else
             write(lun_livei,FMT=3) whea1
             write(lun_livei,FMT=3) whea2//'RESTART ALYA FOR PROBLEM: '//trim(title)
             write(lun_livei,FMT=3) whea1
          end if
          if( mapps > 1 ) then
             write(lun_livei,FMT=3) 'RUNNING ALYA WITH ANOTHER CODE'
             write(lun_livei,FMT=3,ADVANCE='no') 'NAMES OF CODES: '
             do ii = 1,mapps
                if( ii == mapps ) then
                   write(lun_livei,FMT=3) trim(application_names(ii))
                else
                   write(lun_livei,FMT=3) trim(application_names(ii))//', '
                end if
             end do
          else
             continue
          end if
          call iofile_flush_unit(lun_livei)
          return
       else if((itask==2.and.lun_livei==6).or.(itask==-2.and.lun_livei/=6) ) then
          messa = 'READ PROBLEM DATA'
       else if( itask == 3 ) then
          messa = 'MODULE DATA'
          inews = 1
       else if( itask == 4 ) then
          if(kfl_timco/=-1 ) then
             messa = 'START TIME STEP '//trim(intost(ittim+1_ip))
             inews = 1
          else
             return
          end if
       else if( itask == 5 ) then
          if(nblok>1.and.itcou==1 ) then
             messa = 'START BLOCK '//trim(intost(iblok))
             inews = 1
          else
             return
          end if
       else if( itask == 6 ) then
          if(micou(iblok)>1 ) then
             messa = 'START COUPLING ITERATION '//trim(intost(itcou))
             inews = 1
          else
             return
          end if
       else if( itask == 7 ) then
          if(micou(iblok)>1 ) then
             messa = 'END COUPLING ITERATION '
             isect = isect - 3
          else
             return
          end if
       else if( itask == 8 ) then
          if(nblok>1 ) then
             messa = 'END BLOCK'
             isect = isect - 3
          else
             return
          end if
       else if( itask == 9 ) then
          if(kfl_timco/=-1 ) then
             messa = 'END TIME STEP'
             isect = isect - 3
          else
             return
          end if
       else if( itask == 10 ) then
          messa = 'END PROBLEM'
          isect = isect - 3
       else if( itask == 12 ) then
          messa = 'READ MESH DATA'
       else if( itask == 14 ) then
          if(kfl_ptask==0 ) then
             messa = 'COMPUTE GRAPH'
          else
             messa = 'COMPUTE GRAPH'
          end if
       else if( itask == 15 ) then
          messa = 'SOLVE '//trim(namod(inume))
       else if( itask == 16 ) then
          dumml=intost(inume)
          messa = ' ('//trim(dumml)//')'
       else if( itask == 17 ) then
          messa = 'CHECK ELEMENT  ORDERING'
       else if( itask == 21 ) then
          messa = 'CHECK BOUNDARY ORDERING'
       else if( itask == 22 ) then
          messa = 'CHECK ELEMENT  TYPES'
       else if( itask == 23 ) then
          messa = 'CHECK ELEMENT  CONNECTIVITY'
       else if( itask == 24 ) then
          messa = 'CHECK BOUNDARY TYPES'
       else if( itask == 25 ) then
          messa = 'CHECK BOUNDARY CONNECTIVITY'
       else if( itask == 26 ) then
          messa = 'CHECK CONNECTIVITY BOUNDARY/ELEMENT'
       else if( itask == 29 ) then
          messa = 'READ  BOUNDARY CONNECTIVITY AND TYPE'
       else if( itask == 30 ) then
          messa = 'READ  BOUNDARY CONNECTIVITY, TYPES AND ELEMENT CONNECTIVITY'
       else if( itask == 34 ) then
          messa = 'OUTPUT MESH'
       else if( itask == 35 ) then
          messa = 'READ  BOUNDARY CONNECTIVITY'
       else if( itask == 36 ) then
          messa = 'READ  BOUNDARY TYPES'
       else if( itask == 37 ) then
          messa = 'READ  BOUNDARY/ELEMENT CONNECTIVITY'
       else if( itask == 38 ) then
          messa = 'READ  MESH FROM BINARY FILE'
       else if( itask == 39 ) then
          messa = 'WRITE MESH IN BINARY FILE'
       else if( itask == 40 ) then
          messa = 'COMPUTE EXTERIOR NORMALS'
       else if( itask == 41 ) then
          messa = 'COMPUTE MASS MATRIX'
       else if( itask == 42 ) then
          ! THIS IS THE GENERAL LIVINF, TO WRITE ANYTHING YOU NEED ON SCREEN
          messa=trim(message)
       else if( itask == 43 ) then
          messa = 'SOLVE '//trim(namod(inume))
       else if( itask == 44 ) then
          messa = 'SOLVE '//trim(namod(inume))
       else if( itask == 45 ) then
          messa = 'CREATE BOUNDARY SET MESH'
       else if( itask == 46 ) then
          messa = 'CREATE PLANE MESH'
       else if( itask == 47 ) then
          messa=trim(namod(modul))//': OPEN  '//trim(message)//' FILE DYNAMIC COUPLING'
       else if( itask == 48 ) then
          messa=trim(namod(modul))//': CLOSE '//trim(message)//' FILE DYNAMIC COUPLING'
       else if( itask == 49 ) then
          messa = 'CONSTRUCT MESH DATA OF SELECTED EXAMPLE'
       else if( itask == 50 ) then
          messa = 'WARNINGS HAVE BEEN FOUND IN MODULE '//trim(namod(modul))
       else if( itask == 51 ) then
          messa = ''//trim(namod(inume))//': READ DATA'
       else if( itask == 52 ) then
          messa = 'END MODULE DATA'
          isect = isect - 3
       else if( itask == 53 ) then
          messa = ''//trim(namod(inume))//': INITIAL SOLUTION'
       else if( itask == 56 ) then
          messa = 'SOLVE '//trim(namod(inume))
       else if( itask == 58 ) then
          dumml=intost(inume)
          messa = 'CHECK OPENMP. MAX # THREADS= '//trim(dumml)
       else if( itask == 59 ) then
          messa=trim(namod(modul))//': '//trim(message)
       else if( itask == 60 ) then
          messa = 'GENERATE CARTESIAN MESH'
          inews = 1
       else if( itask == 61 ) then
          messa = 'COMPUTE LOCAL NORMAL BASIS'
       else if( itask == 62 ) then
          messa = 'COMPUTE SYMMETRIC GRAPH'
       else if( itask == 63 ) then
          messa = 'END GENERATE CARTESIAN MESH'
          isect = isect - 3
       else if( itask == 64 ) then
          mess1=intost(ioutp(1))
          mess2=intost(ioutp(2))
          mess3=intost(ioutp(3))
       else if( itask == 65 ) then
          messa = 'READ  HANGING NODES FROM BINARY FILE'
       else if( itask == 66 ) then
          messa = 'READ  IMMERSED BOUNDARY FROM BINARY FILE'
       else if( itask == 67 ) then
          if( inume == -2  ) then
             messa = 'KERNEL: READ LAGRANGIAN PARTICLE RESTART FILE'
          else if( modul == -1  ) then
             messa = 'KERNEL: READ IB RESTART FILE'
          else if( modul == 0  ) then
             messa = 'KERNEL: READ RESTART FILE'
          else
             messa = ''//trim(namod(modul))//': READ RESTART FILE'
          end if
       else if( itask == 68 ) then
          if( inume == -2  ) then
             messa = 'KERNEL: WRITE LAGRANGIAN PARTICLE RESTART FILE'
          else if( modul == -1  ) then
             messa = 'KERNEL: WRITE IB RESTART FILE'
          else if( modul == 0  ) then
             messa = 'KERNEL: WRITE RESTART FILE'
          else
             messa=trim(namod(modul))//': WRITE RESTART FILE'
          end if
       else if( itask == 69 ) then
          messa = 'IMMBOU: COMPUTE IMMERSED BOUNDARY PROPERTIES'
       else if( itask == 70 ) then
          messa=trim(namod(modul))//': COMPUTE FORCES AND TORQUES ON IB'
       else if( itask == 71 ) then
          messa = 'IB GAUSS POINTS FALL ON SUBDOMAIN INTERFACE'
       else if( itask == 72 ) then
          messa = 'SOLVE EULER EQUATIONS FOR IB'
       else if( itask == 73 ) then
          messa = 'ADAPT MESH'
       else if( itask == 75 ) then
          messa = 'INITIAL SOLUTION'
          inews = 1
       else if( itask == 76 ) then
          messa = 'END INITIAL SOLUTION'
          isect = isect - 3
       else if( itask == 77 ) then
          messa = ''//trim(namod(modul))
       else if( itask == 78 ) then
          messa = 'SAVE ELEMENT DATA BASE'
       else if( itask == 79 ) then
          messa = ''//trim(namod(inume))//': '//trim(message) ! End of time step
       else if( itask == 80 ) then
          if(kfl_paral<=0)&
               messa = 'RENUMBER NODES'
       else if( itask == 81 ) then
          messa = ''//trim(namod(modul))//': '//trim(message) ! End of time step
       else if( itask == 82 ) then
          messa = 'IMMERSED BOUNDARY GAUSS POINTS'
       else if( itask == 83 ) then
          messa = 'DIVIDE MESH'
       else if( itask == 84 ) then
          messa = 'RECONSTRUCT BOUNDARY/ELEMENT CONNECTIVITY'
       else if( itask == 85 ) then
          messa = 'END ADAPT MESH'
          isect = isect - 3
       else if( itask == 86 ) then
          mess1=intost(nelem)
          mess2=intost(npoin)
          mess3=intost(nboun)
          messa = 'NEW MESH: ELEMENTS= '//trim(mess1)//', NODES= '//trim(mess2)//', BOUNDARIES= '//trim(mess3)
       else if( itask == 88 ) then
          messa = 'WRITE NEW GEOMETRY IN FILES'
       else if( itask == 89 ) then
          messa = 'MARK ELEMENTS TO BE ADAPTED'
       else if( itask == 90 ) then
          messa = 'COMPUTE GEOMETRICAL NORMALS'
       else if( itask == 95 ) then
          messa = 'CHECK PROJECTION'
       else if( itask == 96 ) then
          messa = 'IMMERSED BOUNDARY FRINGE NODES'
       else if( itask == 97 ) then
          continue
       else if( itask == 98 ) then   ! similar to 69 but for ale
          messa = 'ALEFOR: COMPUTE RIGID BODY PROPERTIES'
       else if( itask == 99 ) then   ! similar to 69 but for solidz
          messa = 'SOLIDZ: COMPUTE RIGID BODY PROPERTIES'
       else if( itask == 160 ) then
          messa = ' ('
       else if( itask == 161 ) then
          dumml=intost(inume)
          messa = ' ('//trim(dumml)//','
       else if( itask == 162 ) then
          dumml=intost(inume)
          messa=trim(dumml)//')'
       else if( itask == 163 ) then
          dumml=intost(inume)
          messa=trim(dumml)//','
       else if( itask == 164 ) then
          messa = ')'
       else if( itask == 165 ) then
          messa=trim(message)
       else if( itask == 166 ) then
          messa=trim(message)
       else if( itask == 201 ) then
          isect = isect-3
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//'CRITICAL TIME OVERSHOOT'
          !        write(lun_livei,FMT='(a)') whea2//wspac(1:isect-3)//'RESETTING TIME STEP '//trim(intost(ittim))
          write(lun_livei,FMT='(a,es13.6)') whea2//repeat(' ',max(0_ip,isect))//'GOING BACK TO TIME t= ',cutim
          messa=trim('RESTARTING TIME STEP '//trim(intost(ittim)))
          isect = isect+3
       else if( itask == 202 ) then
          messa = 'CALL INITIALIZE COPROCESSING'
       else if( itask == 203 ) then
          messa = 'CALL COPROCESSING'
       else if( itask == 204 ) then
          messa = 'CALL FINALIZE COPROCESSING'
       else if( itask == 999 ) then
          !write(6,3) ' '
          messa = 'ABORTED IN |---->   '//trim(message)
       else if( itask == 1000 ) then
          !messa = 'FINISHED NORMALLY'
          messa = 'CALCULATIONS CORRECT'
       end if

       !-------------------------------------------------------------------
       !
       ! Generic messages: to be used preferently
       !
       !-------------------------------------------------------------------

       if( itask == 0  ) then
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(message)

       else if( itask == -2  ) then
          !write(lun_livei,FMT='(a,$)') whea2//repeat(' ',max(0_ip,isect))//trim(message)
          write(lun_livei,FMT='(a)',advance='no') whea2//repeat(' ',max(0_ip,isect))//trim(message)

       else if( itask == -3  ) then
          write(lun_livei,FMT='(a)') trim(message)

       else if( itask == -4 ) then
          inews = 1
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(message)

       else if( itask == -5 ) then
          isect = isect - 3
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(message)

       else if( itask == -6 ) then
          inews = 1
          messa = intost(inume)
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '//trim(messa)

       else if( itask == -7 ) then

          write(lun_livei,FMT='(a,i7)') whea2//repeat(' ',max(0_ip,isect))//trim(message)//': ',inume

       else if( itask == -8 ) then

          write(lun_livei,FMT='(a,a1,i7,a1,a,a1,i7)') whea2//repeat(' ',max(0_ip,isect))//trim(coutp(1)),&
               ' ',ioutp(1),' ',trim(coutp(2)),' ',ioutp(2)

       else if( itask == -9 ) then

          if( inume < 0 ) then
             write(lun_livei,FMT='(a,i7)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message),inume
          else if( inume >= 10**6 ) then
             write(lun_livei,FMT='(a,i7)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message),inume
          else if( inume >= 10**5 ) then
             write(lun_livei,FMT='(a,i6)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message),inume
          else if( inume >= 10**4 ) then
             write(lun_livei,FMT='(a,i5)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message),inume
          else if( inume >= 10**3 ) then
             write(lun_livei,FMT='(a,i4)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message),inume
          else if( inume >= 10**2 ) then
             write(lun_livei,FMT='(a,i3)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message),inume
          else if( inume >= 10 ) then
             write(lun_livei,FMT='(a,i2)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message),inume
          else
             write(lun_livei,FMT='(a,i1)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message),inume
          end if

       else if( itask == -10 ) then

          write(lun_livei,FMT='(a,a1,i7,a1,a,a1,i7)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//&
               trim(coutp(1)),' ',ioutp(1),' ',trim(coutp(2)),' ',ioutp(2)

       else if( itask == -11  ) then
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message)

       else if( itask == -12  ) then
          write(lun_livei,FMT='(a)') ' '
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(message)

       else if( itask == -13 ) then
          isect = isect - 3
          messa = intost(inume)
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '//trim(messa)

       else if( itask == -14 ) then
          messa = intost(inume)
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '//trim(messa)

       else if( itask == -15 ) then

          write(lun_livei,FMT='(a,a1,i7)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//&
               trim(coutp(1)),' ',ioutp(1)

       else if( itask == -16 ) then

          write(lun_livei,FMT='(a,e12.6)') whea2//repeat(' ',max(0_ip,isect))//trim(namod(modul))//': '//trim(messaGE)//'= ',routp(1)

       else if( itask == -17  ) then
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//'!!! WARNING: '//trim(message)

       else if( itask == -18  ) then
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//'!!! REPORT: '//trim(message)

       else if( itask == -19 ) then
          inews = 1
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//'START '//trim(message)

       else if( itask == -20 ) then
          isect = isect - 3
          write(lun_livei,FMT='(a)') whea2//repeat(' ',max(0_ip,isect))//'END '//trim(message)

          !-------------------------------------------------------------------
          !
          ! Specific messages
          !
          !-------------------------------------------------------------------

       else if( itask == 4 ) then
          write(lun_livei,FMT='(a)',advance='no') whea2//repeat(' ',max(0_ip,isect))//trim(messa)

       else if( itask == 15.or.itask==44.or.itask==43.or.itask==56 ) then
          write(lun_livei,FMT='(a)',advance='no') whea2//repeat(' ',max(0_ip,isect))//trim(messa)
          !write(lun_livei,FMT='(a)',advance='no') whea2//repeat(' ',max(0_ip,isect))//trim(messa)

       else if( itask == 165.or.itask==161.or.itask==160.or.itask==163 ) then
          write(lun_livei,FMT='(a)',advance='no') trim(messa)
          !write(lun_livei,FMT='(a)',advance='no') trim(messa)

       else if( itask == 164.or.itask==16 ) then
          write(lun_livei,FMT='(a)') trim(messa)

       else if( itask == 166 ) then
          write(lun_livei,FMT=3) whea2//repeat(' ',max(0_ip,isect))//' '

       else if( itask == 17.or.itask==21 ) then
          write(lun_livei,FMT='(a)',advance='no') whea2//repeat(' ',max(0_ip,isect))//trim(messa)

       else if( itask == 18 ) then
          write(lun_livei,FMT='(a,es13.6)') ', t= ',cutim

       else if( itask == 19 ) then
          write(lun_livei,FMT=1,advance='no') inume

       else if( itask == 20 ) then
          write(lun_livei,FMT=*)

       else if( itask == 97 ) then
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
          !mess1 = intost(ioutp(1))//trim(wbyte(1))
          !mess2 = intost(ioutp(2))//trim(wbyte(2))
          !mess3 = intost(ioutp(3))//trim(wbyte(3))
          !messa = ' (NELEM= '//trim(mess1)//', NPOIN= '//trim(mess2)//', NBOUN= '//trim(mess3)//')'
          !write(lun_livei,FMT=3) whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '//adjustl(trim(messa))

          messa = ' (NELEM= '//trim(mess1)//', NPOIN= '//trim(mess2)//', NBOUN= '//trim(mess3)//')'
          write(lun_livei,FMT=5) whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '&
               //' (NELEM= ',(routp(1)),wbyte(1) &
               //', NPOIN= ',(routp(2)),wbyte(2) &
               //', NBOUN= ',(routp(3)),wbyte(3) &
               //')'

          !mess1 = intost(ioutp(1))
          !mess2 = intost(ioutp(2))
          !mess3 = intost(ioutp(3))
          !messa = ' (NELEM= '//trim(mess1)//', NPOIN= '//trim(mess2)//', NBOUN= '//trim(mess3)//')'
          !write(lun_livei,FMT=3) whea2//repeat(' ',max(0_ip,isect))//trim(message)//' '//adjustl(trim(messa))

       else if( itask == 98  ) then
          write(lun_livei,FMT='(a,e16.8E3)') whea2//repeat(' ',max(0_ip,isect))//trim(message),routp(1)

       else if( itask == 99  ) then
          write(lun_livei,FMT='(a,i5)') whea2//repeat(' ',max(0_ip,isect))//trim(message),ioutp(1)

       else if( itask == 999 ) then
          write(lun_livei,FMT=3) whea1
          write(lun_livei,FMT=3) whea2//repeat(' ',max(0_ip,isect))//adjustl(trim(messa))
          write(lun_livei,FMT=3) whea1

       else if( itask == 1000 ) then
          write(lun_livei,FMT=3) whea1
          write(lun_livei,FMT=3) whea2//repeat(' ',max(0_ip,isect))//adjustl(trim(messa))
          write(lun_livei,FMT=3) whea1

       else if( itask == 10000 ) then
          write(lun_livei,FMT=3) whea1
          messa=whea1//' WARNING: '// adjustl(trim(message))
          write(lun_livei,FMT=3) whea2//repeat(' ',max(0_ip,isect))//adjustl(trim(messa))
          write(lun_livei,FMT=3) whea1
          call iofile_flush_unit(lun_livei)

       else if( itask == -1 ) then
          messa = ''//trim(message)
          write(lun_livei,FMT=3) trim(messa)

       else if( itask == 64 ) then
          write(lun_livei,FMT=4) '',trim(mess1),' CELLS MARKED FOR REFINEMENT ON ',&
               &             trim(mess2),' TOTAL CELLS FOR ITER ',trim(mess3)

       else

          write(lun_livei,FMT=3) whea2//repeat(' ',max(0_ip,isect))//trim(messa)

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
