subroutine chm_begmet()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_begmet
  ! NAME 
  !    chm_begmet
  ! DESCRIPTION
  !    This routine prepares a new time step for METEO model
  ! USES
  !    chm_updunk
  ! USED BY
  !    chemic
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none


  if( kfl_meteo_chm < 0) return !Cloud model cannot enter here.

  if( kfl_meteo_chm >= 1 .and. tmete_chm <= oltim ) then
     !
     ! Read meteo file: VELOC_CHM, TEMPE_CHM, DENSI_CHM
     !
     call chm_reamet()
     !
     ! Update terminal velocity each time meteo file is read
     !
     call chm_updvte()
     
  end if
  
  if( lawde_chm == -2 .or. lawte_chm == -2 .or. kfl_advec_chm == -2 ) then
     !
     ! Update terminal velocity 
     !
     call chm_updvte()
  end if
  
  if( kfl_sourc_chm <= -1 .and. tsour_chm <= oltim ) then
     !
     ! Read source file
     !
     call chm_reasou()
  end if
  !
  ! Update boundary conditions
  !
  call chm_updbcs(1_ip)
     

end subroutine chm_begmet
