!-----------------------------------------------------------------------
!> @addtogroup Master
!> @{
!> @file    mod_get_options.f90
!> @author  houzeaux
!> @date    2020-05-11
!> @brief   Get options
!> @details Get Alya options
!-----------------------------------------------------------------------

module mod_get_options

  use def_kintyp_basic, only : ip
  
  implicit none
  private
  
  type getopt_t
     character(1)   :: short  = ''  ! The short option (single character, without the leading dash)
     character(30)  :: long   = ''  ! The long option (without the leading dashes, max 99 characters long)
     integer(4)     :: reqarg = 0_4 ! Argument required? 0-no, 1-yes
     character(50)  :: descr  = ''  ! A (short) description (recommended: <1 screen width; max 999 characters)
  end type getopt_t

  public :: getopt_t
  public :: get_options_help
  public :: get_options_this
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-11
  !> @brief   Help message
  !> @details Help message
  !> 
  !-----------------------------------------------------------------------
  
  subroutine get_options_help(longopts)
    
    type(getopt_t),   intent(in) :: longopts(:)
    character(30)                :: wlong
    character(50)                :: wdescr
    integer(ip)                  :: ii

    write(6,1)
    do ii = 1,size(longopts)
       wlong  = adjustl('--'//adjustl(trim(longopts(ii) % long)))
       wdescr = adjustl(longopts(ii) % descr)
       write(6,'(a,a30,a50)') &
            '    -'//trim(longopts(ii) % short),&
            '  '//adjustl(wlong),&
            wdescr
    end do
    
    write(6,3)
    stop

1   format(&
         & '   ',/,&
         & '   ALYA USAGE:',/,&
         & '   ',/,&
         & '   Alya.x [options] case ',/)
3   format(&
         & '   ',/,&
         & '   ',/,&
         & '   Runs Alya for problem case. ',/,&
         & '   ',/,&
         & '   The following I/O files are located/created in current directory (mod is any activated module extension)',/,&
         & '   * means optional:',/,&
         & '   ',/,&
         & '   (I)    case.dat:                     run data',/,&
         & '   (I)    case.ker.dat:                 kernel data',/,&
         & '   (I)    case.dom.dat:                 mesh data',/,&
         & '   (I*)   case.cou.dat:                 coupling data',/,&
         & '   (I)    case.mod.dat:                 module data',/,&
         & '   ',/,&
         & '   (O)    case.log:                     run log',/,&      
         & '   (O)    case.ker.log:                 kernel log',/,&      
         & '   (O*)   case.mem:                     memory',/,&
         & '   (O*)   case.liv:                     live info',/,&
         & '   (O)    case-partition.par.post.msh   partition mesh in GiD format',/,&
         & '   (O)    case-partition.par.post.res   partition results in GiD format',/,&
         & '   (O)    case-VAR.post.alyabin:        postprocess file of variable VAR',/,&
         & '   (O)    case-VAR.mod.sol              solver information for variable VAR',/,&
         & '   (O)    case-VAR.mod.cso              solver convergence for variable VAR',/,&
         & '   (O)    case.mod.cvg:                 module convergence',/,&      
         & '   (O)    case.mod.rst:                 module restart',/,&  
         & '   ')

  end subroutine get_options_help

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-11
  !> @brief   Get short from long option
  !> @details Long to short 
  !> 
  !-----------------------------------------------------------------------
  
  type(getopt_t) function get_options_this(arg,longopts)

    type(getopt_t),   intent(in) :: longopts(:)   
    character(LEN=*), intent(in) :: arg
    integer(ip)                  :: ii,my_len
    
    get_options_this % short  = ''
    get_options_this % long   = ''
    get_options_this % reqarg = 0_4
    my_len                    = len_trim(arg)
    
    if( len_trim(arg) >= 2 ) then
       if( arg(1:2) == '--' ) then
          do ii = 1,size(longopts)
             if(     trim(longopts(ii) % long ) == trim(arg(3:my_len)) ) then
                get_options_this % short  = longopts(ii) % short
                get_options_this % long   = longopts(ii) % long
                get_options_this % reqarg = longopts(ii) % reqarg
                get_options_this % descr  = longopts(ii) % descr
                return
             end if
          end do          
       else if( arg(1:1) == '-' ) then
          do ii = 1,size(longopts)
             if(     trim(longopts(ii) % short) == trim(arg(2:my_len)) ) then
                get_options_this % short  = longopts(ii) % short
                get_options_this % long   = longopts(ii) % long
                get_options_this % reqarg = longopts(ii) % reqarg
                get_options_this % descr  = longopts(ii) % descr
                return
             end if
          end do
       end if
    end if
 
  end function get_options_this
  
end module mod_get_options
!> @}
