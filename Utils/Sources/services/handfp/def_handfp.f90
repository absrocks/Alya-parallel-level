module def_handfp
!-----------------------------------------------------------------------
!
! A module to handle floating point exceptions
!
!-----------------------------------------------------------------------  

! ia32 winnt compaq visual fortran
  USE DFLIB
  INTERFACE
     FUNCTION handfp (sigid, except)
     !DEC$ ATTRIBUTES C :: handfp
        INTEGER(4) handfp
        INTEGER(2) sigid, except
     END FUNCTION
  END INTERFACE
  integer(4) :: iret

end module def_handfp
