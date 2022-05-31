subroutine tem_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_sendat
  ! NAME
  !    tem_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_temper
  use def_inpout
  use mod_memchk
  use def_kermod,           only : gasco
  use mod_opebcs
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,ki,kr,ir,ii,kfl_ptask_old,dummi
  integer(4)              :: istat

  select case (order)

  case(1_ip)    

     !------------------------------------------------------------------- 
     !
     ! Exchange data read in tem_reaphy, tem_reanut and tem_reaous
     !
     !------------------------------------------------------------------- 

     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
      

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of tem_reaphy variables 
        !
        call iexcha(kfl_timei_tem)
        call iexcha(kfl_advec_tem)
        call iexcha(kfl_joule_tem)
        call iexcha(kfl_radia_tem)
        call iexcha(kfl_sourc_tem)
        call iexcha(kfl_sonum_tem)
        call iexcha(kfl_adiab_tem)
        call iexcha(kfl_condu_tem)
        call iexcha(kfl_exint_tem)
        call iexcha(kfl_inter_tem)
        call iexcha(kfl_dynco_tem)
        call iexcha(kfl_regim_tem)
        call iexcha(kfl_lookg_tem)
        call iexcha(kfl_diven_tem)
        call iexcha(kfl_parti_tem)
        call iexcha(idtem_tem)
        call iexcha(kfl_flux_tem)
        call iexcha(kfl_skews_tem)
        call iexcha(kfl_entropy_tem)
        call iexcha(kfl_code_tem)

        call rexcha(turbu_tem)
        call rexcha(prthe_tem)
        call rexcha(cfi_hmax_tem)
        call rexcha(cfi_hmin_tem)
        call rexcha(prtur_tem)
        call rexcha(scond_tem)
        call rexcha(react_tem)
        call rexcha(sourc_tem)
        call rexcha(gasco)

        do ii = 1,nmate
           call iexcha(lsour_material_tem(ii))
        end do
        do ii = 1,nmate
           do ji = 1,msour_material_tem
              call rexcha(xsour_material_tem(ji,ii))
           end do
        end do
        !
        ! Exchange of tem_reanut variables 
        !        
        call iexcha(kfl_dttyp_tem)
        call iexcha(kfl_ellen_tem)
        call iexcha(kfl_sgsti_tem)
        call iexcha(kfl_sgsno_tem)
        call iexcha(kfl_taust_tem)
        call iexcha(kfl_ortho_tem)
        call iexcha(kfl_limit_tem)
        call iexcha(kfl_shock_tem)
        call iexcha(kfl_tiacc_tem)
        call iexcha(kfl_tibub_tem)
        call iexcha(kfl_assem_tem)
        call iexcha(kfl_posit_tem)
        call iexcha(kfl_negat_tem)
        call iexcha(neule_tem)
        call iexcha(kfl_tisch_tem)
        call iexcha(kfl_normc_tem)
        call iexcha(miinn_tem)
        call iexcha(misgs_tem)
        call iexcha(kfl_sgsac_tem)
        call iexcha(kfl_meshi_tem)
        call iexcha(kfl_discr_tem)
        call iexcha(kfl_sgsli_tem)
        call iexcha(kfl_plepp_tem)
        call iexcha(kfl_explicit_tem)
        call iexcha(kfl_rhs_scal_tem)

        call rexcha(staco_tem(1))
        call rexcha(staco_tem(2))
        call rexcha(staco_tem(3))
        call rexcha(shock_tem)
        call rexcha(safet_tem)
        call rexcha(source_safet_tem)
        call rexcha(sstol_tem)
        call rexcha(cotol_tem)
        call rexcha(relax_tem)
        call rexcha(bemol_tem)
        call rexcha(relsg_tem)
        call rexcha(tosgs_tem) 
        solve_sol => solve(1:)
        call soldef(1_ip)
        !
        ! Exchange data read in tem_reaous
        !
        call posdef(1_ip,dummi)
        call iexcha(kfl_splot_tem)
        call iexcha(kfl_psmat_tem)
        call iexcha(kfl_exacs_tem)
        do ji=1,10
           call iexcha(kfl_viewf_tem(ji))
        end do
        call iexcha(npp_bound_tem)
        call rexcha(avtim_tem)
        do ji=1,nexap_tem
           call rexcha(expar_tem(ji))
        end do
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','tem_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','tem_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call par_broadc()
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call par_broadc()
     
     call memchk(two,istat,mem_modul(1:2,modul),'parin','tem_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','tem_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','tem_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','tem_sendat',0_ip)     
     !
     ! Allocatable arrays: LMATE_TEM and LMATN_TEM
     !
      
   

     kfl_ptask = kfl_ptask_old
     !
     ! Physical properties
     !
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of tem_reaphy variables whose dimensions depend
        ! on what is read in tem_reaphy
        !  
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','tem_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','tem_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call par_broadc()
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'parin','tem_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','tem_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','tem_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','tem_sendat',0_ip)     

      
     
     !------------------------------------------------------------------- 
     !
     ! Variables read in reabcs
     !
     !------------------------------------------------------------------- 

     call spnbcs(tncod_tem)
     !call spnbcs(tgcod_tem)
     call spbbcs(tbcod_tem)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of tem_reabcs variables 
        !
        call iexcha(kfl_conbc_tem) 
        call iexcha(kfl_inidi_tem)
        call iexcha(kfl_inico_tem)
        call iexcha(kfl_intbc_tem)
        call iexcha(npnat_tem)
        call rexcha(delta_tem)   
        call rexcha(initial_tem)
        do jr=1,10
           call iexcha(kfl_funty_tem(jr))
        enddo
        do jr=1,10
           do ir=1,6
              call rexcha(funpa_tem(ir,jr))
           enddo
        enddo
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','tem_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','tem_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call par_broadc()
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'parin','tem_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','tem_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','tem_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','tem_sendat',0_ip)     

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine tem_sendat
