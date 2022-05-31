#ifdef CMAKE

subroutine infmak()
  !------------------------------------------------------------------------
  !****f* info/infmak
  ! NAME
  !    infmak
  ! DESCRIPTION
  !    Output makefile info
  ! USES
  !    outfor
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_master
  implicit none

  coutp(1)  = 'COMPILED BY:'
  coutp(2)  = 'DATE:'
  coutp(3)  = 'HOSTNAME:'
  coutp(4)  = 'ARCHITECTURE:'
  coutp(5)  = 'OPERATING SYSTEM:'
  coutp(6)  = 'END'
  coutp(10) = 'F90:'
  coutp(11) = 'FPPFLAGS:'
  coutp(12) = 'FCFLAGS:'
  coutp(13) = 'FOPT:'
  coutp(14) = 'CSALYA:'
  coutp(15) =  'END'

end subroutine infmak

#endif
