subroutine tur_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_sendat
  ! NAME
  !    tur_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_turbul
  use def_inpout
  use mod_memchk
  use mod_opebcs
  use def_kermod, only     :  cmu_st, kfl_kxmod_ker, kfl_logva
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,ki,dummi
  integer(ip)             :: kfl_ptask_old 
  integer(4)              :: istat

  select case (order)

  case(1_ip)     
     !
     ! Exchange data read in tur_reaphy, tur_reanut and tur_reaous
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
      

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of tur_reaphy variables 
        !
        call iexcha(kfl_model_tur)
        call iexcha(kfl_timei_tur)
        call iexcha(kfl_cotem_tur)
        call iexcha(kfl_advec_tur)
        call iexcha(kfl_colev_tur)
        call iexcha(kfl_ddesm_tur)
        call iexcha(kfl_sasim_tur)
        call iexcha(kfl_rotat_tur)
        call iexcha(kfl_inifi_tur(1))           ! Initial fields 
        call iexcha(kfl_inifi_tur(2))           ! Initial fields 
        call iexcha(kfl_inifi_tur(3))           ! Initial fields
        do jr=1,3
           call iexcha(nfiel_tur(jr))           ! Fields assignement
        end do

        call iexcha(inits_tur)
        call iexcha(nturb_tur)
        do ji=1,nipar_tur
           call iexcha(ipara_tur(ji))
        end do
        call iexcha(lawde_tur)
        call iexcha(lawvi_tur)
        call iexcha(kfl_kxmod_tur)
        call iexcha(kfl_logva_tur)
        call iexcha(kfl_discd_tur)
        call rexcha(boube_tur)
        call rexcha(grnor_tur)
        do jr=1,3
           call rexcha(gravi_tur(jr))
        end do
        call rexcha(prtur_tur)
        do jr=1,nrpar_tur
           call rexcha(param_tur(jr))
        end do
        do jr=1,ncoef_tur
           call rexcha(densi_tur(jr))
        end do
        do jr=1,ncoef_tur
           call rexcha(visco_tur(jr))
        end do
        call rexcha(densa_tur)
        call rexcha(visca_tur)
        call rexcha(cddes_tur)       
        call rexcha(inv_l_max)       
        
        !
        ! Exchange of tur_reanut variables 
        !        
        call iexcha(kfl_weigh_tur)
        call iexcha(kfl_repro_tur)
        call iexcha(kfl_taust_tur)
        call iexcha(kfl_shock_tur)
        call iexcha(kfl_algor_tur)
        call iexcha(kfl_clipp_tur)
        call iexcha(kfl_ellen_tur)
        call iexcha(kfl_relax_tur)
        call iexcha(kfl_tiacc_tur)
        call iexcha(kfl_tisch_tur)
        call iexcha(kfl_normc_tur)
        call iexcha(kfl_walgo_tur)
        call iexcha(kfl_assem_tur)
        call iexcha(kfl_ortho_tur)
        call iexcha(kfl_limit_tur)
        call iexcha(kfl_produ_tur)
        call iexcha(kfl_meshi_tur)

        call iexcha(miinn_tur)
        call iexcha(niter_tur)
        call iexcha(neule_tur)

        call iexcha(kfl_sgsti_tur)
        call iexcha(kfl_sgsno_tur)
        call iexcha(kfl_tibub_tur)
        call iexcha(kfl_sgsac_tur)

        call rexcha(staco_tur(1))
        call rexcha(staco_tur(2))
        call rexcha(staco_tur(3))
        call rexcha(shock_tur)
        call rexcha(sstol_tur)
        call rexcha(safet_tur)
        call rexcha(cotol_tur)
        call rexcha(relax_tur)
        call rexcha(safex_tur)
        call rexcha(safma_tur)
        call rexcha(safeo_tur)
        call rexcha(saflo_tur)
        call rexcha(bemol_tur)

        call rexcha(clipfac_tur)
        call iexcha(kfl_lmaxi_tur)
        solve_sol => solve(1:5)
        call soldef(1_ip)
        !
        ! Exchange data read in tur_reaous
        !
        call posdef(1_ip,dummi)
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','tur_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','tur_sendat',parre)
           if( ISLAVE .or.kfl_ptask==2 ) call par_broadc()
        end if
     end do
     if( IMASTER .and.kfl_ptask/=2 ) call par_broadc()
     
     call memchk(two,istat,mem_modul(1:2,modul),'parin','tur_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','tur_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','tur_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','tur_sendat',0_ip)
     

     kfl_kxmod_ker = kfl_kxmod_tur  ! copy variable to kernel in all processes
     cmu_st      = param_tur(6) ! Overwrites !!! copy variable to kernel variable in all processes
     if (kfl_logva==1) kfl_logva_tur=1 ! copy variable to kernel variable
     if (kfl_logva_tur==1) kfl_logva=1 ! copy variable from kernel variable
     !
     ! Material wake dissipation term
     !
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of tur_reaphy variables 
        !
        if (kfl_discd_tur==1) then
           do ji=1, nmate        
              call iexcha(ldiss_material_tur(ji))
           end do
        end if
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','tur_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','tur_sendat',parre)
           if( ISLAVE .or.kfl_ptask==2 ) call par_broadc()
        end if
     end do
     if( IMASTER .and.kfl_ptask/=2 ) call par_broadc()
     
     call memchk(two,istat,mem_modul(1:2,modul),'parin','tur_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','tur_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','tur_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','tur_sendat',0_ip)
     kfl_ptask = kfl_ptask_old
   case(2_ip)     
     !
     ! Exchange data read in tur_reabcs
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
      
     !
     ! Boundary codes
     !
     call spnbcs(tncod_tur)
     call spgbcs(tgcod_tur)
     call spbbcs(tbcod_tur)
     !
     ! Variables read in tur_reabcs
     !
      
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of tur_reabcs variables 
        !
        call iexcha(kfl_inidi_tur)
        call iexcha(kfl_wallw_tur)
        call iexcha(kfl_infl1_tur)
        call iexcha(kfl_infl2_tur)
        call iexcha(kfl_usrbc_tur)
        call iexcha(kfl_initi_tur)
        do ki = 1,4
           call iexcha(kfl_valbc_tur(ki))
        end do
        call rexcha(delta_tur)
        call rexcha(turin_tur)
        call rexcha(hdiam_tur)
        call rexcha(turle_tur) 
        call rexcha(nutnu_tur)
        call rexcha(rebcs_tur)
        do ki = 1,4
           call rexcha(xinit_tur(ki))
        end do
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','tur_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','tur_sendat',parre)
           if( ISLAVE .or. kfl_ptask==2 ) call par_broadc()
        end if
     end do 

     if( IMASTER .and.kfl_ptask/=2) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'parin','tur_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','tur_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','tur_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','tur_sendat',0_ip)     

  case(3_ip)   !! all this lines should be errased they are no longer called
     !
     ! Arrays computed in tur_addarr
     !
!     do parii=1,2 
!        npari=0
!        nparr=0
!        nparc=0
!        !
!        ! Exchange of tur_reabcs variables 
!        !
!        call iexcha(kfl_walld_tur)
!        call iexcha(kfl_lwnei_tur)
!        call iexcha(kfl_ustar_tur)
!        call iexcha(kfl_grve2_tur)
!        call iexcha(kfl_grsqk_tur)
!        call iexcha(kfl_grono_tur)
!        call iexcha(kfl_greps_tur)
!!        call iexcha(kfl_grk12_tur)
!        call iexcha(kfl_grphi_tur)
!        call iexcha(kfl_vorti_tur)
!        call iexcha(kfl_avvel_tur)
!        call iexcha(kfl_adapt_tur)
!        call iexcha(kfl_fixn6_tur)
!        call iexcha(kfl_fixn8_tur)
        !
        ! Allocate memory for the first pass
        !
!        if(parii==1) then
!           allocate(parin(npari),stat=istat)
!           call memchk(zero,istat,mem_modul(1:2,modul),'parin','tur_sendat',parin)
!           allocate(parre(nparr),stat=istat)
!           call memchk(zero,istat,mem_modul(1:2,modul),'parre','tur_sendat',parre)
!           if( ISLAVE .or. kfl_ptask==2 ) call par_broadc()
!        end if
!     end do

!     if( IMASTER .and. kfl_ptask/=2 ) call par_broadc()

!     call memchk(two,istat,mem_modul(1:2,modul),'parin','tur_sendat',parin)
!     deallocate(parin,stat=istat)
!     if(istat/=0) call memerr(two,'parin','tur_sendat',0_ip)
!     call memchk(two,istat,mem_modul(1:2,modul),'parre','tur_sendat',parre)
!     deallocate(parre,stat=istat)
!     if(istat/=0) call memerr(two,'parre','tur_sendat',0_ip)     

!     call tur_memarr( 5_ip) ! USTAR_TUR: Velocity gradient
!     call tur_memarr( 8_ip) ! KFL_GRK12_TUR, KFL_GRONO_TUR, KFL_GREPS_TUR
!     call tur_memarr( 9_ip) ! GRVE2_TUR: Velocity 2nd order gradients
!     call tur_memarr(10_ip) ! GRSQK_TUR: grad(sqrt(k))
!     call tur_memarr(11_ip) ! GRPHI_TUR: grad(phi)

!     if( INOTMASTER.or.kfl_ptask/=2 ) then
        !
        ! WALLD_TUR
        !
!        if(kfl_walld_tur/=0.and.kfl_walgo_tur==0) then
!           if( INOTMASTER ) call tur_memarr(1_ip)
!           strre =  'WALLD_TUR'
!           call vocabu(NPOIN_REAL_1DIM,0_ip,0_ip)
!           parr1 => walld_tur
!           call par_mygather()
!           if( IMASTER ) call tur_memarr(7_ip)
!        end if  

!     end if

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine tur_sendat
