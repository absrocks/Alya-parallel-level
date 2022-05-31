subroutine nsi_sendat(order)
  !-----------------------------------------------------------------------
  !****f* nastin/nsi_sendat
  ! NAME
  !    nsi_sendat
  ! DESCRIPTION
  !    This routine exchange NASTIN data 
  ! USES
  ! USED BY
  !    nsi_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_nastin
  use def_inpout
  use mod_memchk
  use mod_opebcs
  use def_kermod, only       :  gasco

  use mod_vortex, only: vortex_sendat   

  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,ii,dummi
  integer(4)              :: istat

  select case (order)

  case(1_ip)     
     !
     ! Exchange data read in nsi_reaphy, nsi_reanut and nsi_reaous
     !
      
     strre='nsi_reaphy_nsi_reanut_nsi_reaous'
     strin='nsi_reaphy_nsi_reanut_nsi_reaous'
     strch='nsi_reaphy_nsi_reanut_nsi_reaous'
     do parii = 1,2 
        npari = 0
        nparr = 0
        nparc = 0

        !----------------------------------------------------------------
        !
        ! Exchange of nsi_reaphy variables 
        !
        !----------------------------------------------------------------

        call iexcha(kfl_timei_nsi)
        call iexcha(kfl_advec_nsi)
        call iexcha(kfl_convection_type_nsi)       
        call iexcha(kfl_fvfua_nsi)
        call iexcha(kfl_fvful_nsi)
        call iexcha(kfl_cotem_nsi)
        call iexcha(kfl_grtur_nsi)
        call iexcha(kfl_visco_nsi)
        call iexcha(kfl_colev_nsi)
        call iexcha(kfl_regim_nsi)
        call iexcha(kfl_dynco_nsi)
        call iexcha(kfl_prthe_nsi)
        call iexcha(kfl_surte_nsi)
        call iexcha(kfl_force_nsi)
        call iexcha(kfl_mfrco_nsi)
        call iexcha(kfl_bnods_nsi)
        call iexcha(kfl_hydro_gravity_nsi)
        call iexcha(kfl_anipo_nsi)
        call iexcha(kfl_fscon_nsi)
        call iexcha(kfl_hydro_interface_nsi)
        call iexcha(mfrse_nsi)
        call iexcha(nbnod_nsi)
        call iexcha(nbval_nsi)
        call iexcha(nbtim_nsi)
        do jr=1,2
           call iexcha(nfiel_nsi(jr))           ! Fields assignement
        end do

        call rexcha(nbtdt_nsi)
        call rexcha(fcons_nsi)
        call rexcha(fvins_nsi)
        call rexcha(grnor_nsi)
        call rexcha(fvnoa_nsi)
        call rexcha(fanoa_nsi)
        call rexcha(fvnol_nsi)
        call rexcha(fanol_nsi)
        call rexcha(bougr_nsi)
        call rexcha(boube_nsi)
        call rexcha(boutr_nsi)
        call rexcha(lowtr_nsi)
        do ji=1,2
           call rexcha(turbu_nsi(ji))
        end do
        call rexcha(heihy_nsi)
        call rexcha(surte_nsi)
        call rexcha(tmass_nsi)
        call rexcha(mfrub_nsi)
        call rexcha(ubpre_nsi)
        call rexcha(mfccf_nsi)

        do ji=1,3
           call rexcha(gravi_nsi(ji))
           call rexcha(gravb_nsi(ji))
           call rexcha(fvdia_nsi(ji))
           call rexcha(fvela_nsi(ji))
           call rexcha(fadia_nsi(ji))
           call rexcha(facca_nsi(ji))
           call rexcha(fvdil_nsi(ji))
           call rexcha(fvell_nsi(ji))
           call rexcha(fadil_nsi(ji))
           call rexcha(faccl_nsi(ji))
           call rexcha(frotc_nsi(ji))
        end do
        call rexcha(centr_nsi)
        do ji=1,6
           call rexcha(fvpaa_nsi(ji))
        end do
        do ji=1,6
           call rexcha(fvpal_nsi(ji))
        end do

        !----------------------------------------------------------------
        !
        ! Exchange of nsi_reanut variables 
        !
        !----------------------------------------------------------------

        call iexcha(kfl_penal_nsi)
        call iexcha(kfl_prepe_nsi)
        call iexcha(kfl_dttyp_nsi)
        call iexcha(kfl_ellen_nsi)
        call iexcha(kfl_relax_nsi)
        call iexcha(kfl_relap_nsi)
        call iexcha(kfl_sgsco_nsi)
        call iexcha(kfl_sgsti_nsi)
        call iexcha(kfl_sgsac_nsi)
        call iexcha(kfl_sgsli_nsi)
        call iexcha(kfl_sgscp_nsi)
        call iexcha(kfl_shock_nsi)
        call iexcha(kfl_tiacc_nsi)
        call iexcha(kfl_normc_nsi)
        call iexcha(kfl_refer_nsi)
        call iexcha(kfl_linea_nsi)
        call iexcha(kfl_tisch_nsi)
        call iexcha(kfl_algor_nsi)
        call iexcha(kfl_predi_nsi)
        call iexcha(kfl_taush_nsi)
        call iexcha(kfl_ellsh_nsi)
        call iexcha(kfl_updpr_nsi)
        call iexcha(kfl_intpr_nsi)
        call iexcha(kfl_assem_nsi)
        call iexcha(kfl_asbou_nsi)
        call iexcha(kfl_taust_nsi)
        call iexcha(kfl_stabi_nsi)
        call iexcha(kfl_limit_nsi)
        call iexcha(kfl_trres_nsi)
        call iexcha(kfl_prtre_nsi)
        call iexcha(kfl_matdi_nsi)
        call iexcha(kfl_intfo_nsi)
        call iexcha(kfl_press_nsi)
        call iexcha(kfl_bubbl_nsi)
        call iexcha(momod(modul) % miinn)
        call iexcha(kfl_stain_nsi)
        call iexcha(kfl_immer_nsi)
        call iexcha(kfl_grad_div_nsi)
        call iexcha(misgs_nsi)
        call iexcha(npica_nsi)
        call iexcha(itinn(modul))
        call iexcha(neule_nsi)
        call iexcha(kfl_savco_nsi)
        call iexcha(kfl_meshi_nsi)
        call iexcha(kfl_corre_nsi)
        call iexcha(kfl_sosch_nsi)
        call iexcha(kfl_modfi_nsi)
        call iexcha(kfl_expco_nsi)
        call iexcha(kfl_addpr_nsi)
        call iexcha(kfl_grvir_nsi)
        call iexcha(kfl_hydro_nsi)
        call iexcha(kfl_update_hydro_nsi)
        call iexcha(kfl_hydro_interface_nsi)
        call iexcha(mitri_nsi)
        call iexcha(kfl_adres_nsi)
        call iexcha(kfl_incre_nsi)
        call iexcha(kfl_nota1_nsi)
        call iexcha(kfl_ini_ts_guess_order_nsi)
        call iexcha(kfl_vector_nsi)
        call iexcha(kfl_press_stab_nsi)
        call iexcha(kfl_stop_by_wit_nsi)
        call iexcha(kfl_massm_nsi)
        
        call rexcha(penal_nsi)
        call rexcha(prepe_nsi)
        call rexcha(dtcri_nsi)
        do jr=1,4
           call rexcha(staco_nsi(jr))
        end do
        call rexcha(shock_nsi)
        call rexcha(safet_nsi)
        call rexcha(bemol_nsi)
        call rexcha(sstol_nsi)
        call rexcha(cotol_nsi)
        call rexcha(resid_nsi)
        call rexcha(resip_nsi)
        call rexcha(weigh_nsi)
        call rexcha(relax_nsi)
        call rexcha(relap_nsi)
        call rexcha(relsg_nsi)
        call rexcha(tosgs_nsi)
        call rexcha(strec_nsi)
        call rexcha(dampi_nsi)
        call rexcha(epsht_nsi)
        call rexcha(epstr_nsi)
        call rexcha(xfree_nsi)
        call rexcha(safex_nsi)
        call rexcha(adres_nsi)
        call rexcha(toler_nsi)
        call rexcha(safma_nsi)
        call rexcha(safeo_nsi)
        call rexcha(saflo_nsi)
        call rexcha(gamma_nsi)
        call rexcha(fsrot_nsi)

        solve_sol => solve(1:)
        call soldef(1_ip)

        !----------------------------------------------------------------
        !
        ! Exchange data read in nsi_reabcs
        !
        !----------------------------------------------------------------

        call iexcha(kfl_confi_nsi)
        call iexcha(kfl_local_nsi)
        call iexcha(kfl_conbc_nsi)
        call iexcha(kfl_initi_nsi)
        call iexcha(kfl_inico_nsi)
        call iexcha(kfl_inipr_nsi)
        call iexcha(kfl_nopen_nsi)
        call iexcha(kfl_cadan_nsi)
        call iexcha(kfl_syntu_nsi)
        call iexcha(kfl_aiobo_nsi)
        call iexcha(neddy_nsi)
        call iexcha(itebc_nsi)
        call iexcha(nodpr_global_nsi)
        call iexcha(exfpr_nsi)
        call iexcha(kfl_imppr_nsi)
        call iexcha(kfl_divcorrec_nsi)
        do ii = 1,size(kfl_flow_rate_codes_nsi,KIND=ip)
           call iexcha(kfl_flow_rate_codes_nsi(ii))
        end do
        do ii = 1,size(kfl_flow_rate_stfn_nsi,KIND=ip)
           call iexcha(kfl_flow_rate_stfn_nsi(ii))
        end do
        do ii = 1,size(kfl_flow_rate_tfn_nsi,KIND=ip)
           call iexcha(kfl_flow_rate_tfn_nsi(ii))
        end do
        do ii = 1,size(kfl_flow_rate_normal_nsi,KIND=ip)
           call iexcha(kfl_flow_rate_normal_nsi(ii))
        end do
        call rexcha(delta_nsi)
        call rexcha(relbc_nsi)
        call rexcha(valpr_nsi)
        call rexcha(mulpr_nsi)
        call rexcha(hydro_nsi)
        call rexcha(velin_nsi(1))
        call rexcha(velin_nsi(2))
        call rexcha(velin_nsi(3))
        call rexcha(poise_nsi(1))
        call rexcha(poise_nsi(2))
        call rexcha(poise_nsi(3))
        call rexcha(poise_nsi(4))
        call rexcha(poise_nsi(5))
        call rexcha(poise_nsi(6))
        call rexcha(divcorrec_nsi)
        do ii = 1,size(flow_rate_values_nsi,KIND=ip)
           call rexcha(flow_rate_values_nsi(ii))
           do ji=1,ndime
              call rexcha(flow_rate_normals_nsi(ji,ii))
           end do
        end do

        call vortex_sendat()  

        !----------------------------------------------------------------
        !
        ! Exchange data read in nsi_reaous
        !
        !----------------------------------------------------------------

        call posdef(1_ip,dummi)
        call iexcha(kfl_exacs_nsi)
        call iexcha(kfl_exfix_nsi)
        call iexcha(kfl_inert_nsi)
        call iexcha(kfl_psmat_nsi)
        do ji=1,10
           call rexcha(expar_nsi(ji))
        end do
        call rexcha(cloth_nsi)
        call rexcha(metab_nsi)
        call rexcha(wetme_nsi)
        call rexcha(ambie_nsi)
        call rexcha(radia_nsi)
        call rexcha(relat_nsi)
        call rexcha(avtim_nsi)
        call rexcha(entim_nsi)   ! probably not needed
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','nsi_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','nsi_sendat',parre)
           if( ISLAVE ) call par_broadc()
        end if
     end do

     if( IMASTER ) call par_broadc()
     call memchk(two,istat,mem_modul(1:2,modul),'parin','nsi_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','nsi_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','nsi_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','nsi_sendat',0_ip)
     !
     ! Allocate properties memory for slaves
     !
      
     !
     ! Material force term
     !
     if( kfl_force_nsi == 1 .and. ISLAVE ) call nsi_memphy(5_ip) ! master was allocated before
     if( kfl_bnods_nsi == 1 .and. ISLAVE ) then
        call nsi_memphy(6_ip)
        call nsi_memphy(7_ip)
     end if
     
     do parii = 1,2 
        npari = 0
        nparr = 0
        nparc = 0

        !----------------------------------------------------------------
        !
        ! Exchange properties and material force read in nsi_reaphy
        !
        !----------------------------------------------------------------

        call rexcha(gasco)
        call rexcha(sphea_nsi)
        call rexcha(prthe_nsi)
        !
        ! Material force
        !
        if( kfl_force_nsi == 1 ) then
           do ii = 1,nmate
              call iexcha(lforc_material_nsi(ii))
           end do
           do ii = 1,nmate
              call iexcha(ntabl_nsi(ii))
           end do
           do ii = 1,nmate
              call iexcha(ntabr_nsi(ii))
           end do
           do ii = 1,nmate
              do ji = 1,mforc_material_nsi
                 call rexcha(xforc_material_nsi(ji,ii))
              end do
           end do
           do ii =1, nmate
              do ji =1, mtabl_nsi
                 call rexcha(velta_nsi(ji,ii))
              end do
           end do
           do ii =1, nmate
              do ji =1, mtabl_nsi
                 call rexcha(thrta_nsi(ji,ii))
              end do
           end do
           do ii =1, nmate
              do ji =1, mtabl_nsi
                 call rexcha(powta_nsi(ji,ii))
              end do
           end do
           do ii =1, nmate
              do ji =1, mtabl_nsi
                 call rexcha(veave_nsi(ji,ii))
              end do
           end do
           
           do ii =1, nmate
              do ji =1, mtabl_nsi
                 call rexcha(radiu_nsi(ji,ii))
              end do
           end do
           do ii =1, nmate
              do ji =1, mtabl_nsi
                 call rexcha(forcn_nsi(ji,ii))
              end do
           end do
           do ii =1, nmate
              do ji =1, mtabl_nsi
                 call rexcha(forct_nsi(ji,ii))
              end do
           end do
           
        end if
        !
        ! exchange boundary nodes list
        !
        if( kfl_bnods_nsi == 1 ) then
          do ii =1,nbnod_nsi
             call rexcha(bntab_nsi(ii,1))
             call rexcha(bntab_nsi(ii,2))
             call rexcha(bntab_nsi(ii,3))
          end do
          do ii =1,nbval_nsi*nbtim_nsi
             call rexcha(bnval_nsi(ii,1))
             call rexcha(bnval_nsi(ii,2))
             call rexcha(bnval_nsi(ii,3))
          end do
        end if          
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','nsi_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','nsi_sendat',parre)
           if( ISLAVE ) call par_broadc()
        end if
     end do

     if( IMASTER ) call par_broadc()
     call memchk(two,istat,mem_modul(1:2,modul),'parin','nsi_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','nsi_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','nsi_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','nsi_sendat',0_ip)   

     call spnbcs(tncod_nsi)
     call spgbcs(tgcod_nsi)
     call spbbcs(tbcod_nsi)

  end select

  npari = 0
  nparr = 0
  nparc = 0

end subroutine nsi_sendat
