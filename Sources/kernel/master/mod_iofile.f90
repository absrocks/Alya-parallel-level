!-----------------------------------------------------------------------
!> @defgroup IO_Toolbox
!> Toolbox for IO, like open and close files
!> @{
!> @file    mod_iofile.f90
!> @author  houzeaux
!> @date    2018-04-11
!> @brief   Do things with file
!> @details Do some operations with units and files
!>
!-----------------------------------------------------------------------

module mod_iofile

  use def_kintyp_basic, only : ip,rp,lg
  use def_master,       only : file_opened
  use def_master,       only : kfl_reawr
  use def_master,       only : lun_outpu
  use def_master,       only : intost

  implicit none
  private
  
  interface iofile_flush_unit
     module procedure iofile_flush_unit_4,&
          &           iofile_flush_unit_8
  end interface iofile_flush_unit

  public :: iofile                  ! Do plenty of things with lots of arguments!
  public :: iofile_open_unit        ! Open a unit
  public :: iofile_close_unit       ! Close a unit
  public :: iofile_available_unit   ! Look for an available unit
  public :: iofile_flush_unit       ! Flush a unit
  public :: iofile_opened           ! Inquire if unit is opened
  public :: iofile_append_tag       ! Append a tag to a file (e.g. MPI rank)
  public :: iofile_number_of_lines  ! Computes the number of lines in a file until blank line is found
  public :: iofile_file_exists      ! If a file exists
  public :: iofile_restart_run      ! Define options for a restart run
  public :: iofile_normal_run       ! Define options for a normal run
  public :: iofile_create_directory ! Create a directory
  public :: iofile_delete_file      ! Delete a file if it exists

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Delete a file
  !> @details Delete a file if it exists
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_delete_file(wfile,nunit)

    character(len=*),           intent(in) :: wfile
    integer(ip),      optional, intent(in) :: nunit
    integer(4)                             :: nunit4
    
    if( iofile_file_exists(wfile) ) then
       if( present(nunit) ) then
          nunit4 = int(nunit,4)
       else
          nunit4 = int(iofile_available_unit(10000_ip),4)
       end if
       open(UNIT=nunit4,file=trim(wfile),status='old')
       close(UNIT=nunit4,status='DELETE')
    end if
    
  end subroutine iofile_delete_file
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Inquire if file exists
  !> @details Inquire if file exists
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function iofile_file_exists(wfile)

    character(len=*), intent(in) :: wfile

    inquire(FILE=trim(wfile),EXIST=iofile_file_exists)

  end function iofile_file_exists

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Inquire if unit is opened
  !> @details Inquire if unit is opened
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function iofile_opened(lunit)

    integer(ip), intent(in) :: lunit
    integer(4)              :: unit4

    unit4=int(lunit,4)
    inquire(unit=unit4,opened=iofile_opened)

  end function iofile_opened

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Flush a unit
  !> @details Flush a unit
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_flush_unit_ip(nunit)

    integer(ip), intent(in) :: nunit
    integer(4)              :: iostat4,nunit4
    integer(ip)             :: iostat

    nunit4 = int(nunit,4)
    flush(UNIT=nunit4,IOSTAT=iostat4,ERR=100)
    return

    iostat = int(iostat4,ip)
100 call runend('IOFILE_FLUSH_UNIT: COULD NOT FLUSH UNIT '//intost(nunit)//' WITH ERROR= '//intost(iostat))
    
  end subroutine iofile_flush_unit_ip
  
  subroutine iofile_flush_unit_4(nunit)

    integer(4),  intent(in) :: nunit
    integer(ip)             :: nunit_ip

    nunit_ip = int(nunit,ip)
    call iofile_flush_unit_ip(nunit_ip)

  end subroutine iofile_flush_unit_4
 
  subroutine iofile_flush_unit_8(nunit)

    integer(8),  intent(in) :: nunit
    integer(ip)             :: nunit_ip

    nunit_ip = int(nunit,ip)
    call iofile_flush_unit_ip(nunit_ip)

  end subroutine iofile_flush_unit_8
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Available unit
  !> @details Look for an available unit
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function iofile_available_unit(START_AT)

    integer(ip), optional, intent(in) :: START_AT
    integer(4)                        :: ioerr,unit4,munit4
    integer(4)                        :: unit4_start
    logical(lg)                       :: opened

    munit4 = 10000_4

    if( present(START_AT) ) then
       unit4_start = int(START_AT,4)
    else
       unit4_start = 90_4
    end if
    
    do unit4 = unit4_start,munit4
       inquire(unit=unit4,opened=opened,iostat=ioerr)
       if( ioerr /= 0 )  cycle
       if( .not. opened ) exit
    end do

    if( unit4 > munit4 ) then
       iofile_available_unit = 0_ip
    else
       iofile_available_unit = int(unit4,ip)
    end if

  end function iofile_available_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Open unit
  !> @details Open a file
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_open_unit(nunit,files,messa,stato,formo,posio,conver,crash_if_cannot_open,IOSTAT)

    integer(ip),      intent(in)              :: nunit 
    character(len=*), intent(in)              :: files
    character(len=*), intent(in),    optional :: messa
    character(len=*), intent(in),    optional :: stato
    character(len=*), intent(in),    optional :: formo
    character(len=*), intent(in),    optional :: posio
    character(len=*), intent(in),    optional :: conver
    logical(lg),      intent(in),    optional :: crash_if_cannot_open
    integer(ip),      intent(inout), optional :: IOSTAT
    logical(lg)                               :: crash

    if( present(IOSTAT) ) then
       crash = .false.
    else
       crash = .true.
    end if
    if( present(crash_if_cannot_open) ) crash = crash_if_cannot_open

    if( present(IOSTAT) ) then
       if( crash ) then
          call iofile(0_ip,nunit,files,messa,stato,formo,posio,conver,IOSTAT=IOSTAT)
       else
          call iofile(7_ip,nunit,files,messa,stato,formo,posio,conver,IOSTAT=IOSTAT)
       end if
    else
       if( crash ) then
          call iofile(0_ip,nunit,files,messa,stato,formo,posio,conver)
       else
          call iofile(7_ip,nunit,files,messa,stato,formo,posio,conver)
       end if
    end if

  end subroutine iofile_open_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-11
  !> @brief   Close a file unit
  !> @details Close a file unit
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_close_unit(nunit,files,messa,IOSTAT,STATUS)

    integer(ip),      intent(in)              :: nunit 
    character*(*),    intent(in),    optional :: files
    character*(*),    intent(in),    optional :: messa
    integer(ip),      intent(inout), optional :: IOSTAT
    character(len=*), intent(in),    optional :: STATUS
    integer(4)                                :: iostat4
    integer(4)                                :: nunit4
    integer(4)                                :: lun_outpu4
    
    nunit4 = int(nunit,4)
    
    if(present(STATUS)) then
       close(UNIT=nunit4,status=trim(STATUS),IOSTAT=iostat4,ERR=100)
    else
       close(UNIT=nunit4,IOSTAT=iostat4,ERR=100)
    end if

    if( iostat4 /= 0 ) then
       if( trim(messa) == 'RESTART' ) then
          lun_outpu4 = int(lun_outpu,4)
          if( present(files) ) then
             write(lun_outpu4,101) 'COULD NOT OPEN THEN FOLLOWING RESTART FILE: '//trim(files)
          else
             write(lun_outpu4,101) 'COULD NOT OPEN THEN FOLLOWING RESTART FILE (UNKNOWN NAME)'
          end if
       else
          if( present(files) .and. present(messa) ) then
             call runend('ERROR WHEN CLOSING THE '//trim(messa)//' FILE: '//adjustl(trim(files)))
          else
             call runend('ERROR WHEN CLOSING FILE')
          end if
       end if
    end if

    return
    
100 continue
    if( present(IOSTAT) ) then
       IOSTAT = int(iostat4,ip)
    else
       call runend('IOFILE_FLUSH_UNIT: COULD NOT CLOSE UNIT '//intost(nunit)//' WITH ERROR= '//intost(int(iostat4,ip)))
    end if
    
101 format(5x,'WARNING: ',a)
    
  end subroutine iofile_close_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-11
  !> @brief   Main routine of this module
  !> @details Do all the operations with units and files
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile(itask,nunit,files,messa,stato,formo,posio,conver,lgexi,IOSTAT)

    integer(ip),      intent(in)            :: itask  !< Task
    integer(ip),      intent(in),  optional :: nunit  !< Unit
    character(len=*), intent(in),  optional :: files  !< File name
    character(len=*), intent(in),  optional :: messa  !< Message for warnings and errors
    character(len=*),              optional :: stato  !< Status
    character(len=*),              optional :: formo  !< Form
    character(len=*),              optional :: posio  !< Position
    character(len=*),              optional :: conver !< Convert
    logical(lg),                   optional :: lgexi
    integer(ip),      intent(out), optional :: IOSTAT
    integer(4)                              :: ioerr4
    integer(4)                              :: nunit4
    character(7)                            :: statu
    character(6)                            :: posit
    character(11)                           :: forma
    integer(4)                              :: lun_outpu4

    nunit4=int(nunit,4)

    select case (itask)

    case ( 0_ip )
       !
       ! Open unit
       !
       if( present(files) ) then
          if(present(stato)) then
             statu=stato
          else
             statu='unknown'
          end if
          if(present(formo)) then
             forma=formo
          else
             forma='formatted'
          end if
          if(present(posio)) then
             posit=posio
          else
             posit='asis'
          end if

#ifdef BIG_ENDIAN
          !
          ! Forces all to big endian
          !
          open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,position=posit,convert='BIG_ENDIAN')
#else
          open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,position=posit)
#endif   
          !
          ! Error when opening the file
          !
          if(ioerr4/=0) then
             file_opened = .false.
             if( present(messa) ) then
                call runend('ERROR WHEN OPENING THE '//trim(messa)//' FILE: '//adjustl(trim(files)))
             else
                call runend('ERROR WHEN OPENING FILE: '//adjustl(trim(files)))                
             end if
          else
             file_opened = .true.
          end if
       else
          call runend('IOFILE: FILE NAME NOT GIVEN!')
       end if

    case ( 2_ip )
       !
       ! Close unit
       !
       if(present(stato)) then
          statu=stato
          close(nunit4,status=statu,iostat=ioerr4)
       else
          close(nunit4,iostat=ioerr4)
       end if

       if(ioerr4/=0) then
          if(trim(messa)=='RESTART') then
             lun_outpu4 = int(lun_outpu,4)
             write(lun_outpu4,101) 'COULD NOT OPEN THEN FOLLOWING RESTART FILE: '//trim(files)
          else
             call runend('ERROR WHEN CLOSING THE '//trim(messa)//' FILE: '//adjustl(trim(files)))
          end if
       end if

    case ( 3_ip )
       !
       ! Delete file
       !
       close(nunit4,status='delete',iostat=ioerr4)
       if(ioerr4/=0) then
          call runend('ERROR WHEN CLOSING THE '//trim(messa)//' FILE: '//adjustl(trim(files)))
       end if

    case ( 4_ip )
       !
       ! Check if file exist
       !
       if(present(formo)) then
          forma=formo
       else
          forma='formatted'
       end if
       if(present(posio)) then
          posit=posio
       else
          posit='asis'
       end if
       open(nunit4,file=adjustl(trim(files)),status='old',form=forma,iostat=ioerr4,position=posit)

       if( ioerr4 /= 0 ) then
          file_opened = .false.
          kfl_reawr   = -abs(kfl_reawr)
       else
          file_opened = .true.
          close(nunit4)
       end if

    case ( 7_ip )
       !
       ! Open unit but do not crash if cannot
       !
       if(present(stato)) then
          statu=stato
       else
          statu='unknown'
       end if
       if(present(formo)) then
          forma=formo
       else
          forma='formatted'
       end if
       if(present(posio)) then
          posit=posio
       else
          posit='asis'
       end if

#ifdef BIG_ENDIAN
       ! forces all to big endian
       open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,position=posit,convert='BIG_ENDIAN')
#else
       if( present(conver) ) then
          open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,position=posit)!,convert=conver)
       else
          open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,position=posit)
       end if
#endif   
       if(ioerr4/=0) then
          file_opened = .false.
       else
          file_opened = .true.
       end if

    end select
    !
    ! Format
    !
    if( present(IOSTAT) ) IOSTAT = int(ioerr4,ip)

101 format(5x,'WARNING: ',a)

  end subroutine iofile

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Append a tag to a file
  !> @details Append a tag to a file (e.g. MPI rank)
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_append_tag(filen,tag,fileo)

    character(len=*), intent(inout)          :: filen
    integer(ip),      intent(in)             :: tag
    character(len=*), intent(out),  optional :: fileo
    integer(ip)                              :: fdot,flen
    
    fdot = index(filen,'.')
    flen = len_trim(filen)
    if( fdot == 0 ) fdot = flen
    if( present(fileo) ) then
       fileo = trim(filen(1:fdot-1)//'-'//trim(intost(tag))//filen(fdot:flen))
    else
       filen = trim(filen(1:fdot-1)//'-'//trim(intost(tag))//filen(fdot:flen))
    end if
    
  end subroutine iofile_append_tag

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Number of lines
  !> @details Gives the number of lines of a file
  !> 
  !-----------------------------------------------------------------------

  function iofile_number_of_lines(nunit) result(nlines)

    integer(ip), intent(in)  :: nunit
    integer(4)               :: nunit4
    integer(4)               :: ierro
    integer(ip)              :: nlines
    character(10)            :: str
    
    nunit4 = int(nunit,4)
    nlines = 0
    rewind(nunit)
    do
       read(nunit,FMT='(A)',iostat=ierro,end=10) str
       if( len_trim(str) == 0 ) goto 10
       if( ierro/= 0 ) exit
       nlines = nlines + 1
    end do
    
10  rewind(nunit)

  end function iofile_number_of_lines
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   File option for a restart run
  !> @details File option for a restart run
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_restart_run(stato,formo,posio,FORMATTED,UNFORMATTED)

    character(len=*), intent(out)          :: stato
    character(len=*), intent(out)          :: formo
    character(len=*), intent(out)          :: posio
    logical(lg),      intent(in), optional :: FORMATTED
    logical(lg),      intent(in), optional :: UNFORMATTED
    logical(lg)                            :: if_formatted

    if_formatted = .true.
    if( present(FORMATTED) )   if_formatted = FORMATTED
    if( present(UNFORMATTED) ) if_formatted = .not. UNFORMATTED

    stato = 'old'
    posio = 'append'

    if( if_formatted ) then
       formo = 'formatted'
    else
       formo = 'unformatted'
    end if
   
  end subroutine iofile_restart_run
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   File option for a normal run
  !> @details File option for a normal run
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_normal_run(stato,formo,posio,FORMATTED,UNFORMATTED)

    character(len=*), intent(out)          :: stato
    character(len=*), intent(out)          :: formo
    character(len=*), intent(out)          :: posio
    logical(lg),      intent(in), optional :: FORMATTED
    logical(lg),      intent(in), optional :: UNFORMATTED
    logical(lg)                            :: if_formatted

    if_formatted = .true.
    if( present(FORMATTED) )   if_formatted = FORMATTED
    if( present(UNFORMATTED) ) if_formatted = .not. UNFORMATTED

    stato = 'unknown'
    posio = 'asis'
    
    if( if_formatted ) then
       formo = 'formatted'
    else
       formo = 'unformatted'
    end if

  end subroutine iofile_normal_run

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Create a directory
  !> @details Create a directory
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_create_directory(dirname)

    character(len=*), intent(inout) :: dirname

    call execute_command_line ('mkdir -p '//trim(dirname))
    
  end subroutine iofile_create_directory
  
end module mod_iofile
!> @}
