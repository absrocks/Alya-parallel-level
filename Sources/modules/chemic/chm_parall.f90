subroutine chm_parall(order)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_parall
  ! NAME
  !    chm_sendat
  ! DESCRIPTION
  !    This routine exchanges data 
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_chemic
  use mod_memchk
  use mod_opebcs
  use def_kermod,         only : gasco
  use def_kermod,         only : lookup_fw
  use mod_interp_tab,     only : tab_par_exchange
  use mod_interp_tab,     only : fw_par_exchange 
  use mod_interp_tab,     only : fw_allocate
  use mod_interp_tab,     only : tab_init_fw
#ifdef CANTERA
  use cantera
#endif
  use mod_chm_sectional_soot_model, only: nsect_ssm,nclas_ssm,nspec_ssm 
  use mod_chm_sectional_soot_model, only: RhoC_ssm,Vmax_ssm,RadSoot_ssm 
  use mod_chm_sectional_soot_model, only: indexS,gasCoupling_ssm,nPAH,nSurf, &
                                          idNucl,idCond,idSurf,              &
                                          ID_NUCL,ID_COND,ID_COAG,ID_SURF

  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ii,ji,ki,ki2,kfl_ptask_old,iclas,jclas
  integer(ip)             :: icoef,dummi,ind
  integer(4)              :: istat
    
  
  if( ISEQUEN ) return
  
  select case (order)

  case(1_ip)    
     !
     ! Exchange data read in chm_reaphy, chm_reanut and chm_reaous
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
      

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0

        !----------------------------------------------------------------
        !
        ! Exchange of chm_reaphy variables 
        !
        !----------------------------------------------------------------

        call iexcha(kfl_model_chm)
        call iexcha(kfl_timei_chm)
        call iexcha(kfl_advec_chm)
        call iexcha(kfl_diffu_chm)
        call iexcha(kfl_transport_chm)
        call iexcha(kfl_norma_chm)
        call iexcha(kfl_cotur_chm)
        call iexcha(kfl_radia_chm)
        call iexcha(kfl_pfa_chm)
        call iexcha(kfl_z_chm)
        call iexcha(kfl_key_chm)
        call iexcha(kfl_freq_chm)
        call iexcha(kfl_premix_chm)
        call iexcha(kfl_varYc_chm)
        call iexcha(kfl_varZ_chm)
        call iexcha(kfl_ufpv_chm)
        call iexcha(kfl_tdac_write_chm)
        call iexcha(kfl_heat_loss_chm)
        call iexcha(kfl_lookg_chm)
        call iexcha(kfl_tab_fw_chm)
        call iexcha(kfl_post_fw_chm)
        call iexcha(kfl_entropy_chm)
        call iexcha(kfl_spec_name_chm)
        call iexcha(kfl_htran)
        call iexcha(kfl_spray_chm)
        call iexcha(kfl_field_chm(1))
        call iexcha(kfl_field_chm(2))
        call iexcha(kfl_field_chm(3))
        !
        ! Soot model
        !
        call iexcha(kfl_soot_chm)

        do ind =1,5
           call iexcha(indexS(ind))
        end do

        call iexcha(nPAH)
        call iexcha(nSurf)

        do ind =1,10
           call iexcha(idNucl(10))
           call iexcha(idCond(10))
           call iexcha(idSurf(10))
        end do

        call iexcha(ID_NUCL)
        call iexcha(ID_COND)
        call iexcha(ID_COAG)
        call iexcha(ID_SURF)
        call iexcha(nsect_ssm)
        call iexcha(gasCoupling_ssm)

        call rexcha(RhoC_ssm)
        call rexcha(Vmax_ssm)
        call rexcha(RadSoot_ssm)

        ! 
        ! End sectional soot model variables 
        ! 

        call iexcha(nclas_chm)
        call iexcha(nspec_chm)
        call iexcha(nspec_ssm)
        call iexcha(nclas_ssm)
        call iexcha(nsect_ssm)
        call iexcha(nsize_mech_name)
        call iexcha(nsize_red)
        call rexcha(prthe_chm)
        if (kfl_coupl(ID_CHEMIC,ID_NASTIN) == 0 ) then
           call rexcha(prthe)
        end if
        call rexcha(radwt_chm)
        call rexcha(dac_crit_chm)
        call rexcha(dac_cor_chm)
        call rexcha(sorad_chm)
        call rexcha(socen_chm(1))
        call rexcha(socen_chm(2))
        call rexcha(socen_chm(3))
        call rexcha(surf_tension_chm)
        call rexcha(bf_fuel_chm)
        call rexcha(bo_oxy_chm)
        
        call iexcha(nreac_chm)

        call iexcha(kfl_droplet_id_chm)
        call iexcha(droplet_postprocess_frequency_chm)
        call rexcha(levelSet_threshold_chm)
        call rexcha(droplet_compactness_limit_chm)
        call rexcha(droplet_max_diameter_chm)
        call rexcha(droplet_h_factor_chm)
        
        call iexcha(stofu_chm(1))
        call iexcha(stofu_chm(2))
        call iexcha(stofu_chm(3))

        call rexcha(gasco)

        nparc = nparc + 5

        call iexcha(kfl_ellen_chm)
        call iexcha(kfl_taust_chm)
        call iexcha(kfl_shock_chm)
        call iexcha(kfl_stabi_chm)
        call iexcha(kfl_limit_chm)
        call iexcha(kfl_tiacc_chm)
        call iexcha(kfl_split_chm)
        call iexcha(kfl_wallc_chm)
        call iexcha(kfl_tibub_chm)
        call iexcha(kfl_tisch_chm)
        call iexcha(kfl_normc_chm)
        call iexcha(kfl_dtcri_chm)
        call iexcha(kfl_negat_chm)
        call iexcha(kfl_posit_chm)
        call iexcha(kfl_warni_chm)
        call iexcha(kfl_temli_chm)

        call rexcha(staco_chm(1))
        call rexcha(staco_chm(2))
        call rexcha(staco_chm(3))
        call rexcha(shock_chm)
        call rexcha(bemol_chm)
        call rexcha(temli_chm)
        call rexcha(cotol_chm)
        call rexcha(safet_chm)
        call rexcha(chemical_time_factor)
        call rexcha(cutof_chm)
        call rexcha(sstol_chm)
        call rexcha(strec_chm)
        call rexcha(dampi_chm)
        call rexcha(epsht_chm)
        call rexcha(epstr_chm)
        call rexcha(dtmin_chm)
        call rexcha(dtmax_chm)
        call rexcha(relax_chm)

        ! Variables for CMC combustion model
        call iexcha(kfl_weigh_in_eq_CMC_chm)
        call iexcha(kfl_split_CFD_CMC)
        call iexcha(kfl_solve_enth_CMC_chm)
        call iexcha(nZ_CMC_chm)
        call iexcha(nvar_CMC_chm)
        call iexcha(nZ_AMC_CMC_chm)
        call iexcha(nS_AMC_CMC_chm)
        call iexcha(index_N2)
        call rexcha(Zs_CMC_chm)
        call rexcha(Smax_AMC_CMC_chm)
        call rexcha(S_threshold)


        solve_sol => solve
        call soldef(1_ip)

        !----------------------------------------------------------------
        !
        ! Exchange data read in chm_reaous
        !
        !---------------------------------------------------------------- 
        call posdef(1_ip,dummi)
        do ki=1,npara_chm
           call iexcha(ipara_chm(ki))
        end do
        do ki=1,npara_chm
           call rexcha(rpara_chm(ki))
        end do

        !----------------------------------------------------------------
        !
        ! Allocate memory for the first pass
        !
        !----------------------------------------------------------------

        if( parii == 1 ) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','chm_sendat',parin)

           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','chm_sendat',parre)
           if( ISLAVE .or. kfl_ptask==2 ) call par_broadc()
        end if
     end do

     if( IMASTER .and. kfl_ptask/=2 ) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'parin','chm_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','chm_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','chm_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','chm_sendat',0_ip)    

 
     !
     ! Exchange tables
     !
     call tab_par_exchange(table_coords,   table_tab   )
     call tab_par_exchange(posttab_coords, posttab_tab )
     if (kfl_model_chm == 1) then
        call tab_par_exchange(yc_scaling_coords, yc_scaling_tab )
        call tab_par_exchange(h_scaling_coords, h_scaling_tab )
        call tab_par_exchange(chi_shape_coords, chi_shape_tab )
     endif

     if ( kfl_tab_fw_chm == 0) then
         if (ISLAVE) then
            allocate(table_fw)
            call tab_init_fw(table_fw)
         endif
         call fw_par_exchange(1_ip,table_fw)
         table_fw % main_table => table_tab
         call fw_par_exchange(2_ip,table_fw)
         if (ISLAVE) call fw_allocate(2_ip,table_fw,0_ip)
         do ii = 1, table_fw % main_table % ndim
             if (table_fw % kfl_scale(ii) == 1) then 
                 select case (table_fw % main_table % coords(ii) % name)
                 case ('CMEAN','C    ')
                    table_fw % scaling(ii) % tab => yc_scaling_tab
                 case ('CHIST')
                    table_fw % scaling(ii) % tab => chi_shape_tab
                 case ('IMEAN','I    ')
                    table_fw % scaling(ii) % tab => h_scaling_tab
                 end select
             endif
             if (ISLAVE) call fw_allocate(2_ip,table_fw,ii)
         enddo
     else if(kfl_tab_fw_chm > 0) then
         table_fw => lookup_fw(kfl_tab_fw_chm)
     endif

     if ( kfl_post_fw_chm == 0) then
         if (ISLAVE) then
            allocate(posttable_fw)
            call tab_init_fw(posttable_fw)
         endif
         call fw_par_exchange(1_ip,posttable_fw)
         posttable_fw % main_table => posttab_tab
         call fw_par_exchange(2_ip,posttable_fw)
         if (ISLAVE) call fw_allocate(2_ip,posttable_fw,0_ip)
         do ii = 1, posttable_fw % main_table % ndim
             if (posttable_fw % kfl_scale(ii) == 1) then 
                 select case (posttable_fw % main_table % coords(ii) % name)
                 case ('CMEAN','C    ')
                    posttable_fw % scaling(ii) % tab => yc_scaling_tab
                 case ('CHIST')
                    posttable_fw % scaling(ii) % tab => chi_shape_tab
                 case ('IMEAN','I    ')
                    posttable_fw % scaling(ii) % tab => h_scaling_tab
                 end select
             endif
             if (ISLAVE) call fw_allocate(2_ip,posttable_fw,ii)
         enddo
     else if(kfl_post_fw_chm > 0) then
         posttable_fw => lookup_fw(kfl_post_fw_chm)
     endif
     !
     ! Allocatable arrays 
     !
     if( ISLAVE ) then
        call chm_memphy(1_ip) ! Allocate: Class dependent properties
        if (kfl_model_chm == 4) then
           ! Slaves allocate memory for the matrices whose size was obtained by the master
           call chm_memphy(3_ip)
           call chm_memphy(5_ip)
           call chm_memphy(6_ip)
           call chm_memphy(7_ip)
        end if
     end if

     kfl_ptask = kfl_ptask_old
     !
     ! Broadcast LREAC_CHM
     !
     if (kfl_model_chm == 3 .or. kfl_model_chm == 4) then
        if( INOTMASTER ) &
            allocate(character(len=nsize_mech_name) :: mechanism_path)
     endif

     if (kfl_model_chm == 3) then
        if( INOTMASTER ) &
            allocate(character(len=nsize_red) :: Red_spec)
     endif

     !
     ! Allocate Field index in slaves 
     !
     if (kfl_model_chm == 3 .and. kfl_spec_name_chm > 0_ip) then
        if( INOTMASTER ) &
            allocate(Field_ind_chm(kfl_spec_name_chm))
     endif
     !
     ! Physical properties
     !

     do parii=1,2 

        npari=0
        nparr=0
        nparc=0

        !----------------------------------------------------------------
        !
        ! Exchange of chm_reaphy variables whose dimensions depend
        ! on what is read in chm_reaphy
        !
        !----------------------------------------------------------------
        
        !
        ! Finite-rate chemistry or CMC models
        !
        if (kfl_model_chm == 3 .or. kfl_model_chm == 4) then
           call cexcha(nsize_mech_name,mechanism_path)
#ifdef CANTERA
           if( (INOTMASTER) .and. (parii==2)) &
                  gas_chm = importPhase(mechanism_path)
#endif
        end if


        !
        ! Finite-rate chemistry model
        !
        if (kfl_model_chm == 3) then

           call cexcha(nsize_red,Red_spec)

           do ki=1,kfl_spec_name_chm
              call iexcha(Field_ind_chm(ki))
           end do
        end if


        !
        ! Variables from CMC combustion model
        !
        if (kfl_model_chm == 4) then
           do ki = 1,nZ_CMC_chm
               call rexcha(Z_CMC_chm(ki))
               call rexcha(Xnormalized_prof_CMC_chm(ki))
            end do
            do ki = 1,nZ_CMC_chm-1
               call rexcha(diff_Z_CMC_chm(ki))
            end do
            do ki = 1,nZ_CMC_chm
               call rexcha(temp_inert_CMC_chm(ki))
               call rexcha(temp_equil_CMC_chm(ki))
               do ki2 = 1,nclas_chm+1
                  call rexcha(rscal_inert_CMC_chm(ki,ki2))
                  call rexcha(rscal_equil_CMC_chm(ki,ki2))
               end do
            end do

            ! Variables for AMC model
            do ki = 1,nZ_AMC_CMC_chm
               call rexcha(Z_AMC_CMC_chm(ki))
               do ki2 = 1,nS_AMC_CMC_chm
                  call rexcha(Xintegrated_table_AMC_CMC_chm(ki,ki2))
               end do
            end do
            do ki2 = 1,nS_AMC_CMC_chm
               call rexcha(S_AMC_CMC_chm(ki2))
            end do
        end if

        do iclas=1,nspec_chm
           call rexcha(diffu_chm(1,iclas))
           call rexcha(diffu_chm(2,iclas))
           call rexcha(Le_k(iclas))
        end do
        
        do iclas=1,nspec_chm ! Transfer of species type
           call rexcha(speci(iclas)%visco(1))
           call rexcha(speci(iclas)%visco(2))
           call iexcha(speci(iclas)%lawvi)
           call rexcha(speci(iclas)%weigh)
           call rexcha(speci(iclas)%densi(1))
           call rexcha(speci(iclas)%densi(2))
           call rexcha(speci(iclas)%densi(3))
           call rexcha(speci(iclas)%densi(4))
           call rexcha(speci(iclas)%entha(1))
           call rexcha(speci(iclas)%entha(2))
           call rexcha(speci(iclas)%prand)
           call rexcha(speci(iclas)%lewis)
           do ji = 1,5
              call rexcha(speci(iclas)%trang(ji))
           enddo      
           do ji = 1,4     
              do icoef = 1,10
                 call rexcha(speci(iclas)%cpcoe(icoef,ji))
              enddo
           enddo
           do icoef = 1,2
              call rexcha(speci(iclas)%activ(icoef))
           enddo
        end do

        if( parii == 1 ) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','chm_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','chm_sendat',parre)
           if( ISLAVE .or. kfl_ptask==2 ) call par_broadc()
        end if
     end do

     if( IMASTER .and. kfl_ptask/=2 ) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'parin','chm_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','chm_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','chm_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','chm_sendat',0_ip)    

  case(2_ip)
     !-------------------------------------------------------------------
     !
     ! Boundary conditions
     !
     !-------------------------------------------------------------------

     call spnbcs(tncod_chm,INCLUDE_CHARACTER=.false.)
     call spbbcs(tbcod_chm)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        call iexcha(kfl_allcl_chm)
        call iexcha(kfl_fields_scale_chm)
        call rexcha(bo_oxy_chm)
        call rexcha(bf_fuel_chm)
  
        do iclas = 1,nclas_chm
           call iexcha(kfl_initi_chm(iclas))
        end do

        do iclas = 1,nclas_chm
           call rexcha(xinit_chm(iclas,1))
        end do
        do iclas = 1,nclas_chm
           call rexcha(xinit_chm(iclas,2))
        end do

        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','chm_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','chm_sendat',parre)
           if( ISLAVE .or. kfl_ptask==2 ) call par_broadc()
        end if
     end do


     if( IMASTER .and. kfl_ptask/=2 ) call par_broadc()


     call memchk(two,istat,mem_modul(1:2,modul),'parin','chm_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','chm_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','chm_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','chm_sendat',0_ip)     

      

  end select

end subroutine chm_parall
