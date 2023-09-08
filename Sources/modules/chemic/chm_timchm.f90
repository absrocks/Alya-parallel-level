subroutine  chm_timchm(pgaus,gpmas,gpden,gpcon,chtim)  
  !-----------------------------------------------------------------------
  ! NAME 
  !    chm_timchm
  ! DESCRIPTION
  !    Compute chemical time scale
  ! USES
  ! USED BY
  !    chm_updtsc
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_chemic, only      :  cutof_chm,safet_chm,chemical_time_factor

  implicit none
  integer(ip),  intent(in)  :: pgaus
  real(rp),     intent(in)  :: gpmas(pgaus)
  real(rp),     intent(in)  :: gpden(pgaus)
  real(rp),     intent(in)  :: gpcon(pgaus)
  real(rp),     intent(out) :: chtim

  integer(ip)               :: igaus
  
  ! Compute the minimum among the gauss points of
  !
  !   ( rho * Yk ) / ( Chemical Reaction Rate * safety_factor )
  !
  ! This gives the critical time step of the chemical reaction 
  !
  chtim=1.0e-6
  do igaus = 1,pgaus
     if (gpden(igaus)*gpcon(igaus) > cutof_chm) then ! Cutoff concentrations below some level, they are too small to track
        chtim = max(chtim,abs(gpmas(igaus)/(gpden(igaus)*gpcon(igaus)*chemical_time_factor))) !  inverse time
     endif
  enddo
  chtim=1.0_rp/chtim

end subroutine chm_timchm
