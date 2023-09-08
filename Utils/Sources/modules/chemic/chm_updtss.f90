subroutine chm_updtss()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updtss
  ! NAME 
  !    chm_updtss
  ! DESCRIPTION
  !    This routine computes the time step size
  ! USED BY
  !    chm_timste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ADR,            only : FROM_CRITICAL
  implicit none 
  real(rp) :: dtmin

  if( ADR_chm(1) % kfl_time_integration /= 0 ) then

     dtmin = huge(1.0_rp)

     if( kfl_dttyp_chm == 1 ) then
        !
        ! Time step based on critical time step
        !
        if (kfl_model_chm == 5_ip) then
           call chm_updtcc_cfi(dtmin)
        else
           call chm_updtcc(dtmin)
        endif        

     else if( kfl_dttyp_chm == 2 ) then
        !
        ! Adaptive time step
        !
        call chm_updtsf(dtmin)

     end if

     dtcri_chm = dtmin
     if( dtcri_chm /= 0.0_rp ) dtinv_chm = 1.0_rp/(dtcri_chm*safet_chm)
     if( kfl_timco == 1 )     then
        dtinv = max(dtinv,dtinv_chm)
     endif
  end if

end subroutine chm_updtss
