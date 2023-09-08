function Handfp (signum, excnum)
!-----------------------------------------------------------------------
!
! Floating point exception handler 
!
!-----------------------------------------------------------------------
  implicit none
  integer(2) :: signum, excnum, Handfp
  
  Handfp= 1

  call runend('HANDFP: CANNOT USE THIS FLOATING POINT EXCEPTION HANDLER: COMPILE THE SERVICE')

end function Handfp
