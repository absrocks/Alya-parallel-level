subroutine nsa_chkpar
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_chkpar
  ! NAME 
  !    nsa_chkpar
  ! DESCRIPTION
  !    This routine checks the master-slave-alone status.
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_nastal
  implicit none


  !checking...
  if (kfl_paral > 0) then           ! slaves
     imaster= .false.
     islave = .true.
     iloner = .false.
  else if (kfl_paral == 0) then     ! master
     imaster= .true.
     islave = .false.
     iloner = .false.
  else if (kfl_paral < 0) then      ! alone
     imaster= .false.
     islave = .false.
     iloner = .true.     
  end if
  
  weparal= .false.
  if (imaster .or. islave) weparal= .true.

end subroutine nsa_chkpar
