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

    !!DMM-ADR if(kfl_timco/=2) then
    !!DMM-ADR    dtinv_chm=dtinv
    !!DMM-ADR    if(momod(modul) % kfl_stead==1) dtinv_chm = 0.0_rp
    !!DMM-ADR    if(kfl_timei_chm==0) dtinv_chm = 0.0_rp  

    !!DMM-ADR    if(kfl_tisch_chm==1) then
    !!DMM-ADR       !
    !!DMM-ADR       ! Trapezoidal rule: Euler iterations
    !!DMM-ADR       !
    !!DMM-ADR       if(ittim<=neule_chm) then
    !!DMM-ADR          kfl_tiacc_chm=1
    !!DMM-ADR       else
    !!DMM-ADR          kfl_tiacc_chm=kfl_tiaor_chm
    !!DMM-ADR       end if
    !!DMM-ADR       if(kfl_tiacc_chm==2) dtinv_chm = 2.0_rp*dtinv_chm
    !!DMM-ADR    else
    !!DMM-ADR       !
    !!DMM-ADR       ! BDF scheme: increase integration order at each time step
    !!DMM-ADR       !
    !!DMM-ADR       kfl_tiacc_chm=min(kfl_tiaor_chm,ittim)
    !!DMM-ADR       call parbdf(kfl_tiacc_chm,pabdf_chm)
    !!DMM-ADR    end if
    !!DMM-ADR end if

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
