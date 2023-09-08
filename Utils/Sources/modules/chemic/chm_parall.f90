subroutine chm_parall(order)
  !-----------------------------------------------------------------------
  !****f* partis/chm_parall
  ! NAME
  !    chm_sendat
  ! DESCRIPTION
  !    This routine exchange data 
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
  use def_kermod, only     : gasco
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,ki,kfl_ptask_old,iclas,jclas
  integer(ip)             :: ireac,icoef,dummi
  integer(4)              :: istat

  if( ISEQUEN ) return

  select case (order)

  case(1_ip)     
     !
     ! Exchange data read in chm_reaphy, chm_reanut and chm_reaous
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     call Parall(29_ip)
     if(ISLAVE) call chm_memphy(6_ip) ! Allocate: CFI table

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
        call iexcha(kfl_stagg_chm)
        call iexcha(kfl_timei_chm)
        call iexcha(kfl_advec_chm)
        call iexcha(kfl_diffu_chm)
        call iexcha(kfl_react_chm)
        call iexcha(kfl_lhsas_chm)
        call iexcha(kfl_corve_chm)
        call iexcha(kfl_norma_chm)
        call iexcha(kfl_sourc_chm)
        call iexcha(kfl_meteo_chm)
        call iexcha(kfl_tfles_chm)
        call iexcha(kfl_cotur_chm)
        call iexcha(kfl_radia_chm)
        call iexcha(kfl_tucfi_chm)
        call iexcha(kfl_uncfi_chm)
        call iexcha(kfl_field_chm(1))
        call iexcha(kfl_field_chm(2))
        call iexcha(kfl_benme_chm)

        call iexcha(lawte_chm)
        call iexcha(lawde_chm)
        call iexcha(lawvt_chm)
        call iexcha(nclas_chm)
        call iexcha(nodes_chm)

        call rexcha(sourc_chm)
        call rexcha(radwt_chm)
        call rexcha(sorad_chm)
        call rexcha(socen_chm(1))
        call rexcha(socen_chm(2))
        call rexcha(socen_chm(3))
        call rexcha(tfles_chm)
        call rexcha(flbet_chm)
        call rexcha(table_cfi(1)%fmima(1))
        call rexcha(table_cfi(1)%fmima(2))
        call rexcha(table_cfi(1)%imima(1))
        call rexcha(table_cfi(1)%imima(2))
        do iclas = 1,15 
          do jclas = 1,2
            call rexcha(table_cfi(1)%inval(jclas,iclas))
          end do
        end do
        call iexcha(nreac_chm)
        call iexcha(ncoef_chm)
        call iexcha(table_cfi(1)%nvcfi(1))
        call iexcha(table_cfi(1)%nvcfi(2))
        call iexcha(table_cfi(1)%nvcfi(3))
        call iexcha(table_cfi(1)%nvcfi(4))
        call iexcha(table_cfi(1)%nvcfi(5))
        call iexcha(table_cfi(1)%ndcfi)
        call iexcha(table_cfi(1)%nrcfi)
        call iexcha(table_cfi(1)%nccfi)
        call iexcha(table_cfi(1)%nfcfi)
        call iexcha(table_cfi(1)%nclas)
        call iexcha(stofu_chm(1))
        call iexcha(stofu_chm(2))
        call iexcha(stofu_chm(3))
        call iexcha(kfl_arreh_chm)

        call rexcha(denma_chm)
        call rexcha(radbi_chm)
        call rexcha(gasco)

        do ki = 1,npart_chm
           call rexcha(temma_chm(ki))
        end do
        call rexcha(boltz_chm)
        call rexcha(strat_chm)

        call iexcha(sponge_chm)
        call rexcha(visco_factor_chm)
        call rexcha(visco_axis(1))
        call rexcha(visco_axis(2))
        call rexcha(visco_axis(3))
        call rexcha(visco_range(1))
        call rexcha(visco_range(2))

        nparc = nparc + 5
        if( parii == 2 .and. IMASTER )  parch(1:5)     = wprob_chm(1:5)
        if( parii == 2 .and. ISLAVE  )  wprob_chm(1:5) = parch(1:5)

        !----------------------------------------------------------------
        !
        ! Exchange of chm_reanut variables 
        ! 
        !----------------------------------------------------------------

        call iexcha(kfl_dttyp_chm)
        call iexcha(kfl_ellen_chm)
        call iexcha(kfl_taust_chm)
        call iexcha(kfl_shock_chm)
        call iexcha(kfl_stabi_chm)
        call iexcha(kfl_limit_chm)
        call iexcha(kfl_tiacc_chm)
        call iexcha(kfl_assem_chm)
        call iexcha(kfl_wallc_chm)
        call iexcha(kfl_tibub_chm)
        call iexcha(miinn_chm)
        call iexcha(neule_chm)
        call iexcha(kfl_tisch_chm)
        call iexcha(kfl_normc_chm)
        call iexcha(kfl_coupl_chm)
        call iexcha(kfl_dtcri_chm)
        call iexcha(kfl_dttar_chm)
        call iexcha(kfl_sgsti_chm)
        call iexcha(kfl_negat_chm)
        call iexcha(kfl_posit_chm)
        call iexcha(kfl_warni_chm)
        call iexcha(kfl_meshi_chm)
        call iexcha(kfl_temli_chm)

        call iexcha(kfl_gauss_chm)
        call iexcha(kfl_spite_chm)
        call iexcha(initial_fraction_step_chm)
        call iexcha(max_fixed_point_iterations_chm)

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
        call rexcha(fixed_point_tolerance_chm)
        call rexcha(odeint_tolerance_chm)
        call rexcha(timestep_min_chm)

        solve_sol => solve
        call soldef(1_ip)

        !----------------------------------------------------------------
        !
        ! Exchange data read in chm_reaous
        !
        !----------------------------------------------------------------
        
        call posdef(1_ip,dummi)
        call iexcha(kfl_sized_chm)
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
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','chm_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','chm_sendat',parre)
           if( ISLAVE .or. kfl_ptask==2 ) call Parall(2_ip)
        end if
     end do

     if( IMASTER .and. kfl_ptask/=2 ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','chm_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','chm_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','chm_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','chm_sendat',0_ip)     
     !
     ! Allocatable arrays 
     !
     call Parall(30_ip)
     if( ISLAVE ) then
        call chm_memphy(1_ip) ! Allocate: Class dependent properties
        call chm_memphy(3_ip) ! Allocate: Reactions
        call chm_memphy(5_ip) ! Allocate: CFI table
     end if

     kfl_ptask = kfl_ptask_old
     !
     ! Broadcast LREAC_CHM
     !
     pard1 =  nreac_chm
     pai1p => lreac_chm
     call Parall(37_ip)
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
        if (kfl_model_chm == 5) then
           do iclas = 1,table_cfi(1)%ndcfi
              do jclas = 1,maxval(table_cfi(1)%nvcfi)
                 call rexcha(table_cfi(1)%ivcfi(iclas,jclas))
              end do
           end do
           do iclas=1,table_cfi(1)%nrcfi
             do jclas = 1,table_cfi(1)%nccfi
               call rexcha(table_cfi(1)%table(iclas,jclas))
             end do
           end do
           if (kfl_uncfi_chm == 1) then
              do iclas=1,table_cfi(1)%nvcfi(3)
                 do jclas = 1,3
                    call rexcha(table_cfi(1)%ymass(iclas,jclas))
                 end do
              end do
           end if 
       end if

        do iclas=1,nclas_chm
           call iexcha(lawdi_chm(1,iclas))
           call iexcha(lawdi_chm(2,iclas))
        end do
        do ireac=1,nreac_chm
           do iclas=1,nclas_chm
              call rexcha(stoic_chm(iclas,ireac,1))
              call rexcha(stoic_chm(iclas,ireac,2))
           end do
        end do
        do iclas=1,nclas_chm
           do jclas=1,nclas_chm
              call rexcha(interaction_chm(iclas,jclas))
              call rexcha(interaction_chm(iclas,jclas))
           end do
        end do

        do iclas=1,nclas_chm
           call rexcha(diffu_chm(1,iclas))
           call rexcha(diffu_chm(2,iclas))
        end do
        do ireac=1,nreac_chm
           call rexcha(radiu_chm(ireac))
        end do
        do ireac=1,nreac_chm
           do icoef=1,ncoef_chm
              call rexcha(react_chm(icoef,ireac))
           end do
        end do
        do ireac=1,nreac_chm
           do iclas=1,nclas_chm
              call rexcha(order_chm(iclas,ireac,1))               
              call rexcha(order_chm(iclas,ireac,2))               
           end do
        enddo
        do iclas=1,nclas_chm
           call rexcha(equil_chm(1,iclas))
           call rexcha(equil_chm(2,iclas))
        end do
        do ireac=1,nreac_chm
           do iclas=1,nclas_chm
              call rexcha(effic_chm(iclas,ireac)) ! reaction efficiency for combustion
           enddo
        enddo
        do iclas=1,nclas_chm ! Transfer of species type
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

        do iclas=1,nclas_chm
           call rexcha(diame_chm(iclas))
           call rexcha(rhopa_chm(iclas))
           call rexcha(shape_chm(iclas))
           call rexcha(spher_chm(iclas))
           call rexcha(fract_chm(iclas))
        end do

        if( parii == 1 ) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','chm_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','chm_sendat',parre)
           if( ISLAVE .or. kfl_ptask==2 ) call Parall(2_ip)
        end if
     end do

     if( IMASTER .and. kfl_ptask/=2 ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','chm_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','chm_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','chm_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','chm_sendat',0_ip)     


  case(2_ip)

     !-------------------------------------------------------------------
     !
     ! Boundary conditions
     !
     !-------------------------------------------------------------------

     call spnbcs(tncod_chm)
     call spbbcs(tbcod_chm)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        call iexcha(kfl_allcl_chm)

        do iclas = 1,nspec_chm
           call iexcha(kfl_initi_chm(iclas))
        end do

        do iclas = 1,nclas_chm
           call iexcha(kfl_usrbc_chm(iclas))
        end do
        do iclas = 1,nclas_chm
           call rexcha(xinit_chm(iclas,1))
        end do
        do iclas = 1,nclas_chm
           call rexcha(xinit_chm(iclas,2))
        end do
        do iclas = 1,nclas_chm
           call rexcha(panat_chm(1,iclas))
           call rexcha(panat_chm(2,iclas))
        end do


        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','chm_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','chm_sendat',parre)
           if( ISLAVE .or. kfl_ptask==2 ) call Parall(2_ip)
        end if
     end do

     if( IMASTER .and. kfl_ptask/=2 ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','chm_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','chm_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','chm_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','chm_sendat',0_ip)     

     call Parall(27_ip)

  end select

end subroutine chm_parall
