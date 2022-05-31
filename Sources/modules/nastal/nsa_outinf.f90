subroutine nsa_outinf
!-----------------------------------------------------------------------
!****f* Nastal/nsa_outinf
! NAME 
!    nsa_outinf
! DESCRIPTION
!    This routine writes on the module output files
! USED BY
!    nsa_turnon
!***
!-----------------------------------------------------------------------
  use      def_master

  use      def_nastal

  implicit none
  character(60) :: equat

  if(kfl_paral<=0) then
  end if

end subroutine nsa_outinf
