subroutine ibm_endste()
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_endste
  ! NAME 
  !    ibm_endste
  ! DESCRIPTION
  !    This routine ends a time step of the particle tracking
  ! USES
  !    ibm_cvgunk
  !    ibm_updunk
  !    ibm_output
  ! USED BY
  !    Ibmper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_immbou
  implicit none  
  !
  ! Compute convergence residual of the time evolution
  !
  if( kfl_stead_ibm == 0 .and. kfl_timei_ibm == 1 ) then
     !call ibm_cvgunk(three)
     call ibm_updunk(5_ip)
  end if
  !
  ! Write restart file
  !
  call ibm_restar(2_ip)
  !
  ! If not steady, go on
  !
  if( kfl_stead_ibm == 0 .and. kfl_timei_ibm == 1 .and. kfl_conve(modul) == 1 ) kfl_gotim = 1

end subroutine ibm_endste
