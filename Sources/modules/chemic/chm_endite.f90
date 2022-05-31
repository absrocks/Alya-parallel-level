subroutine chm_endite(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_endite
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
  use mod_chm_finiteRate,         only : chm_getProp_finiteRate
  use mod_chm_finiteRate,         only : chm_IntegrateSource_finiteRate
  use mod_chm_finiteRate,         only : chm_heatRelease_finiteRate
  use mod_chm_operations_CMC,     only : chm_integrate_chem_source_CMC, &
                                         chm_integrate_flow_var_points_CMC, &
                                         chm_integrate_flow_var_gauss_CMC, &
                                         chm_calc_densi_visco_gauss_CMC, &
                                         chm_rho_visco_nodal_project_CMC, &
                                         chm_limit_Yk_CMC, &
                                         chm_global2local_CMC, &
                                         chm_heatRelease_integral_CMC
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: itask, imixf
  real(rp)    :: dt_split

  select case(itask)
     
  case( ITASK_ENDINN )
     !
     ! Cut off over/under shoots
     !
     call chm_clippi()
     !
     ! Compute residual: ||UNKNO(:)-CONCE(:,ICLAS_CHM,1)||
     !     
     call chm_cvgunk(ITASK_ENDINN)
     !
     ! Update unknown: CONCE(:,ICLAS_CHM,1)=UNKNO
     !
     call chm_updunk(ITASK_ENDINN)

  case(ITASK_ENDITE)
     !
     !  Compute residual of coupling iteration + update unknown
     !  This step is called after the module has converged. It is the final step before leaving the module
     !
     call livinf(16_ip,' ',itinn(modul))
     call chm_cvgunk(ITASK_ENDITE) ! Residual: ||CONCE(:,:,2)-CONCE(:,:,1)||
     call chm_updunk(ITASK_ENDITE) ! Update:   CONCE(:,:,2) = CONCE(:,:,1)

     call parari('SUM',0_ip,1_ip,kfl_under_chm)
     if( kfl_under_chm > 0 .and. INOTSLAVE ) then
        call livinf(-9_ip,'UNDERSHOOTS HAVE BEEN FOUND= ',kfl_under_chm)
     endif

     call parari('SUM',0_ip,1_ip,kfl_overs_chm)
     if( kfl_overs_chm > 0 .and. INOTSLAVE ) then 
        call livinf(-9_ip,'OVERSHOOTS HAVE BEEN FOUND= ',kfl_overs_chm)
     endif


     if (kfl_model_chm == 1) then
       !
       ! Read flamelet table for gas phase
       !
       if (kfl_spray_chm == 0 .or. ( kfl_spray_chm /= 0 .and. kfl_premix_chm == 0)) then

          !
          ! Calculate scalar dissipation rate for UFPV model
          !
          if (kfl_ufpv_chm > 0) then
              call chm_post_scalar_dissipation_rate(47_ip)
          endif

          !
          ! Table lookup on Gauss points OR on nodes
          !
          if (kfl_lookg_chm > 0) then
              call chm_gp_reatab()
          else
              call chm_reatab()
          endif

       end if
       call chm_upwmea(ITASK_ENDITE)                         ! wmean(ipoin,1) ==> wmean(ipoin,2)

     elseif (kfl_model_chm == 3) then

       !
       ! Integrate source terms
       !
       if ( kfl_split_chm > 0_ip  ) then
          dt_split = 1.0_rp / dtinv / dble(kfl_split_chm)
          call chm_IntegrateSource_finiteRate(dt_split)
          call chm_heatRelease_finiteRate
       end if

       !
       ! Calculate transport properties
       !
       call chm_getProp_finiteRate()
       call chm_upwmea(ITASK_ENDITE)                         ! wmean(ipoin,1) ==> wmean(ipoin,2)


     elseif (kfl_model_chm == 4) then

       !
       ! Integrate source terms for CMC
       !
       if ( kfl_split_chm > 0_ip  ) then
          dt_split = 1.0_rp / dtinv / dble(kfl_split_chm)
          call chm_integrate_chem_source_CMC(dt_split)
       end if

       call chm_limit_Yk_CMC
       ! Get density and viscosity at Gauss points at all mixt. frac. levels
       do imixf = 1, nZ_CMC_chm
          call chm_global2local_CMC(imixf)
          call chm_calc_densi_visco_gauss_CMC(imixf)
       end do
       ! Do integrations
       call chm_integrate_flow_var_points_CMC
       call chm_integrate_flow_var_gauss_CMC
       call chm_rho_visco_nodal_project_CMC
       call chm_heatRelease_integral_CMC   ! Compute the integral of heat release in the whole domain

       if (kfl_start_CMC_chm == 1)    kfl_start_CMC_chm = 0

     endif


  case(ITASK_INNITE)
      !
      ! After Runge-Kutta substeps:
      ! * Clip under and overshoots
      ! * CONCE(:,ICLAS_CHM,1)=UNKNO
      !
     if ( kfl_negat_chm == 1 .or. kfl_posit_chm == 1) &
          call chm_clippi
     
      call chm_updunk(ITASK_INNITE)

  end select

end subroutine chm_endite
