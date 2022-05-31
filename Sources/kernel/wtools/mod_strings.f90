!-----------------------------------------------------------------------
!> @addtogroup Tools
!> @{
!> @file    mod_strings.f90
!> @author  houzeaux
!> @date    2020-02-27
!> @brief   Strings
!> @details Operations with strings
!-----------------------------------------------------------------------

module mod_strings

  use def_kintyp_basic, only : ip,rp
  implicit none
  private
  
  character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

  public :: upper_case
  public :: lower_case
  public :: integer_to_string
  public :: real_to_string
  public :: string_to_integer
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Integer to string
  !> @details Convert an integer(ip) to a string
  !> 
  !-----------------------------------------------------------------------

  function integer_to_string(integ) result(intost)

    integer(ip),      intent(in)   :: integ
    integer(4)                     :: integ4
    character(len=:), allocatable  :: intost
    character(20)                  :: intaux

    integ4 = int(integ,4)
    write(intaux,*) integ4
    intost = trim(adjustl(intaux))

  end function integer_to_string

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   String to integer
  !> @details Convert a string to an integer
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function string_to_integer(str,stat) result(integ)

    character(len=*),intent(in)            :: str
    integer(ip),     intent(out), optional :: stat
    integer(4)                             :: stat_loc
    
    read(str,*,iostat=stat_loc) integ
    if( present(stat) ) stat = int(stat_loc,ip)
    
  end function string_to_integer

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-27
  !> @brief   Real to string
  !> @details Convert real(rp) to a string
  !> 
  !-----------------------------------------------------------------------

  function real_to_string(realn,REAL_FORMAT) result(retost)

    real(rp)                               :: realn
    character(len=*), intent(in), optional :: REAL_FORMAT
    character(len=:), allocatable          :: retost
    character(20)                          :: reaux
    integer(ip)                            :: ierr

    if( present(REAL_FORMAT) ) then
       write(reaux,REAL_FORMAT,IOSTAT=ierr) realn
       if( ierr /= 0 ) reaux = '0.0'
    else
       write(reaux,'(e19.12)') realn
    end if
    retost = trim(adjustl(reaux))

  end function real_to_string

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Convert to upper case
  !> @details Convert to upper case
  !> 
  !-----------------------------------------------------------------------

  function upper_case (str) result (string)

    character(*),       intent(In) :: str
    character(len(str))            :: string
    integer(ip)                    :: ic, i

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, len_trim(str)
        ic = index(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

  end function upper_case

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-09-29
  !> @brief   Convert to lower case
  !> @details Convert to lower case
  !> 
  !-----------------------------------------------------------------------

  function lower_case (str) result (string)

    character(*),       intent(In) :: str
    character(len(str))            :: string
    integer(ip)                    :: ic, i

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, len_trim(str)
        ic = index(cap, str(i:i))
        if (ic > 0) string(i:i) = low(ic:ic)
    end do

  end function lower_case
  
end module mod_strings
!> @}
