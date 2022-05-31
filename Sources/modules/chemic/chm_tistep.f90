subroutine chm_tistep()
  !-----------------------------------------------------------------------
  !****f* partis/chm_tistep
  ! NAME 
  !    chm_tittim
  ! DESCRIPTION
  !    This routine sets the time step
  ! USES
  ! USED BY
  !    chm_begite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_chemic
  use mod_ADR,    only : ADR_time_strategy
  use mod_outfor, only : outfor
  implicit none
  integer(ip)              :: iclas

  do iclas = 1, nclas_chm
     call ADR_time_strategy(ittim,dtinv,dtinv_old,ADR_chm(iclas))
  end do 

  routp(1) = dtcri_chm
  routp(2) = 0.0_rp
  routp(3) = 0.0_rp
  ioutp(1) = kfl_timei_chm
  ioutp(2) = momod(modul) % kfl_stead
  call outfor(8_ip,lun_outpu,' ')

end subroutine chm_tistep
