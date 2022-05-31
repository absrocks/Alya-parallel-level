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
  use mod_ADR,                   only : FROM_CRITICAL
  use mod_chm_finiteRate,        only : chm_updtcc_finiteRate
  use mod_chm_operations_CMC,    only : chm_updtcc_CMC

  implicit none 
  integer(ip) :: iclas
  real(rp)    :: dtmin

  if( ADR_chm(1) % kfl_time_integration /= 0 ) then

     dtmin = huge(1.0_rp)

     !
     ! Time step based on critical time step
     !
     if (kfl_model_chm == 1) then
        call chm_updtcc_flamLet(dtmin)

     else if (kfl_model_chm == 3) then
        call chm_updtcc_finiteRate(dtmin)

     else if (kfl_model_chm == 4) then
        call chm_updtcc_CMC(dtmin)
     endif        

     dtcri_chm = dtmin

     if( dtcri_chm /= 0.0_rp ) dtinv_chm = 1.0_rp/(dtcri_chm*safet_chm)
     if( kfl_timco == 1 )     then
        dtinv = max(dtinv,dtinv_chm)
     endif

     do iclas = 1,nclas_chm
        ADR_chm(iclas) % dtinv = dtinv_chm
     end do

  end if

end subroutine chm_updtss
