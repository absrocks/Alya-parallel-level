subroutine wav_outlat
!-----------------------------------------------------------------------
!****f* Wavequ/wav_outlat
! NAME 
!    wav_outlat
! DESCRIPTION
!    This routine writes info on the heat equation in latex format
! USED BY
!    wav_turnon
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_wavequ
  implicit none
  character(300) :: equat
  character(20)  :: lvisc
  integer(ip)    :: ierhs

  if(kfl_latex==0) return

end subroutine wav_outlat
      
