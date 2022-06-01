module unitt
  implicit none
  PUBLIC  :: assert_true, assert_false
contains
  subroutine assert_true(my_result,wname,message)  
    logical, intent(in) :: my_result
    character(*), optional :: wname    
    character(*), optional :: message    
    if( .not. my_result ) then
       if( present(message) .and. present(wname) ) then
          write(*,*) 'TEST '//TRIM(wname)//' HAS FAILED IN TEST '//trim(message)
       else
          call abort()
       end if
      end if
  end subroutine assert_true


  subroutine assert_false(my_result,wname,message)
    logical, intent(in) :: my_result
    character(*), optional :: wname    
    character(*), optional :: message    
      if( my_result ) then
       if( present(message) .and. present(wname) ) then
          write(*,*) 'TEST '//TRIM(wname)//' HAS FAILED IN TEST '//trim(message)
       else
          call abort()
       end if
      end if
  end subroutine assert_false
end module unitt

