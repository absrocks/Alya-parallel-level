subroutine exm_chkpar
  !-----------------------------------------------------------------------
  !****f* exmedi/exm_chkpar
  ! NAME 
  !    exm_chkpar
  ! DESCRIPTION
  !    This routine checks the master-slave-alone status.
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_exmedi
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

end subroutine exm_chkpar
