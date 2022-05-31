subroutine ibm_tistep()
!-----------------------------------------------------------------------
!****f* Immbou/ibm_tistep
! NAME 
!    ibm_tistep
! DESCRIPTION
!    This routine sets the time step. 
! USES
! USED BY
!    ibm_begite
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_immbou
  use mod_outfor, only : outfor
  implicit none

  routp(1) = dtcri_ibm
  ioutp(1) = kfl_timei_ibm
  ioutp(2) = kfl_stead_ibm
  routp(2) = 0.0_rp
  routp(3) = 0.0_rp
  call outfor(8_ip,lun_outpu,' ')

end subroutine ibm_tistep
