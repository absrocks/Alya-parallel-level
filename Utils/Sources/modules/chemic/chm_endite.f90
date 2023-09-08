subroutine chm_endite(itask)
  !-----------------------------------------------------------------------
  !****f* partis/chm_endite
  ! NAME 
  !    chm_endite
  ! DESCRIPTION
  !    This routine checks convergence and performs updates of the
  !    temperature  at:
  !    - itask=1 The end of an internal iteration
  !    - itask=2 The end of the internal loop iteration
  ! USES
  !    chm_cvgunk
  !    chm_updunk
  ! USED BY
  !    chm_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_chemic
  use mod_ker_proper 
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: itask

  select case(itask)

  case(1_ip)

     !
     ! Jacobi: Compute residual of internal iteration + update unknown
     ! This option is called after all species have been updated, but still inside intra-chemic interations
     !
     if( kfl_coupl_chm == 1 ) then
        call chm_cvgunk(1_ip) ! Residual:   ||UNKNO(:)-CONCE(:,:,1)||
        call chm_updunk(3_ip) ! Update:     CONCE(:,:,1)=UNKNO
     end if

     ! update other properties
     if( kfl_model_chm == 4 ) then  ! Combustion
        if (kfl_norma_chm > 0) call chm_updunk(9_ip) ! Fix normalization of species
        if ( kfl_gauss_chm == 1_ip) then ! Update eq coefficients at the end of the iteration only
           call chm_upwmea(6_ip) 
           call ker_updpro() 
           call chm_omegak(1_ip,nspec_chm)
        endif
     endif

  case(2_ip)
     !
     !  Compute residual of coupling iteration + update unknown
     !  This step is called after the module has converged. It is the final step before leaving the module
     !  Fractional step combustion model will evolve the reactions here
     !
     if (kfl_reset /= 1) call livinf(16_ip,' ',itinn(modul))
     call chm_cvgunk(2_ip) ! Residual: ||CONCE(:,:,2)-CONCE(:,:,1)||
     call chm_updunk(4_ip) ! Update:   CONCE(:,:,2) = CONCE(:,:,1)

     call parari('SUM',0_ip,1_ip,kfl_under_chm)
     if( kfl_under_chm > 0 .and. INOTSLAVE ) then
        call livinf(-9_ip,'UNDERSHOOTS HAVE BEEN FOUND= ',kfl_under_chm)
     endif

     call parari('SUM',0_ip,1_ip,kfl_overs_chm)
     if( kfl_overs_chm > 0 .and. INOTSLAVE ) then 
        call livinf(-9_ip,'OVERSHOOTS HAVE BEEN FOUND= ',kfl_overs_chm)
     endif

     if (kfl_model_chm == 4) then
       !
       if (kfl_norma_chm > 0) call chm_updunk(9_ip) ! Normalize
       !
       call chm_omegak(1_ip,nspec_chm)              ! Final update Mass source terms
       !
       call chm_upwmea(6_ip)                        
       call chm_upwmea(4_ip)                        
       !
       call chm_heatso()     
     elseif (kfl_model_chm == 5) then
       !
       call chm_reatab()                             ! Update table properties to start doiter in temper with updated coefficients  
       call chm_nsacfi()                             ! send chemical heat to gauss points as needed by nastal
       call chm_upwmea(4_ip)                         ! wmean(ipoin,1) ==> wmean(ipoin,2)
       !
     endif

  case(3_ip)
     !
     ! Gauss-Seidel step: Compute residual of internal iteration + update unknown
     ! This is called inside the module at each GS step and at each intraspecies iteration
     !
     call chm_updunk(10_ip)! Relax update
     call chm_updunk(8_ip) ! Cut off undershoots
     call chm_cvgunk(1_ip) ! Residual:   ||UNKNO(:)-CONCE(:,ICLAS_CHM,1)||
     call chm_updunk(3_ip) ! Update:     CONCE(:,ICLAS_CHM,1)=UNKNO
     !
     ! Update projection
     !
     call chm_solsgs()
 
  end select

end subroutine chm_endite
