subroutine chm_endste()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_endste
  ! NAME 
  !    chm_endste
  ! DESCRIPTION
  !    This routine ends a time step of the transport equation.
  ! USES
  !    chm_cvgunk
  !    chm_updunk
  !    chm_output
  !    chm_restar
  ! USED BY
  !    Chemic
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_chemic
  implicit none

  !
  ! Compute convergence residual of the time evolution (that is,
  ! || c(n,*,*) - c(n-1,*,*)|| / ||c(n,*,*)||) and update unknowns
  ! c(n-1,*,*) <-- c(n,*,*) 
  !
  if(momod(modul) % kfl_stead==0.and.kfl_timei_chm==1) then
     call chm_cvgunk(three)
     call chm_updunk(five)
  end if

  !
  !  If necessary, computes accumulation at "ground"
  !  Ground (surface of accumulation) is defined by the 
  !  kfl_fixno_chm = 2
  !
  if( kfl_model_chm == 2) then
      call chm_accumu()
  end if

  !
  ! In combustion, update variable for other modules
  !
  if (kfl_model_chm == 4) then
     !
     if (kfl_norma_chm > 0) call chm_updunk(9_ip)   ! Normalize
     !
     call chm_omegak(1_ip,nspec_chm)                ! Final update Mass source terms
     !
     call chm_upwmea(5_ip)                          ! Mean molecular weight
     !
     call chm_heatso()                              ! Heat source term
     !
  elseif (kfl_model_chm == 5) then
     call chm_reatab()
     call chm_upwmea(5_ip)                          ! wmean(ipoin,1) ==> wmean(ipoin,3)
  endif
  !
  ! Compute averaged variables
  !
  call chm_averag()

  !
  ! Write restart file
  !
  call chm_restar(2_ip)

  !
  ! If not steady, go on
  !
  if(momod(modul) % kfl_stead==0) then
     if(kfl_timei_chm==1) kfl_gotim = 1
  end if

end subroutine chm_endste
