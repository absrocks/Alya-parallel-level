subroutine dod_outinf
  !-----------------------------------------------------------------------
  !****f* dodeme/dod_outinf
  ! NAME 
  !    dod_outinf
  ! DESCRIPTION
  !    This routine computes some info about Dodeme
  ! USES
  !    
  ! USED BY
  !    dod_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_dodeme
  use mod_outfor, only : outfor
  implicit none
  if(kfl_rstar/=2) then
     call outfor(27_ip,lun_outpu_dod,' ')
  else
     call outfor(11_ip,lun_outpu_dod,' ')
  end if

end subroutine dod_outinf
