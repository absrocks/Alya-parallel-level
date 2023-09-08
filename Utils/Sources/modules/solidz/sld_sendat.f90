!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_sendat.f90
!> @author  Solidz Team
!> @date    November, 2017-Adds doxygen comments
!> @brief   Exchange data
!> @details It is used by sld_turnon
!> @}
!-----------------------------------------------------------------------
subroutine sld_sendat(order)

  use def_kintyp, only : ip, rp
  use def_master, only : IMASTER, ISLAVE
  use def_master, only : nparc, npari, nparr
  use def_master, only : strin, strre
  use def_master, only : parii, parre, parin
  use def_master, only : kfl_ptask, mem_servi, id_parall, kfl_eccty
  use def_domain, only : ndime
  use mod_memchk
  use mod_opebcs
  use def_solidz

  implicit none

  integer(ip), intent(in) :: order                    !> Data read selection
  integer(ip)             :: ji,jr,ki,ii,idumy,ifunc
  integer(ip)             :: kfl_ptask_old
  integer(4)              :: istat

  select case (order)

  case(1_ip)

     !----------------------------------------------------------------
     !
     ! Exchange variables read in sld_reaphy, sld_reanut and sld_reaous
     !
     !----------------------------------------------------------------

     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     call Parall(29_ip)

     do parii=1,2
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of sld_reaphy variables
        !
        call iexcha(kfl_timei_sld)
        call iexcha(kfl_rigid_sld)
        call iexcha(kfl_tange_sld)
        call iexcha(kfl_strai_sld)
        call iexcha(kfl_sdvar_sld)
        call iexcha(kfl_fiber_sld)
        call iexcha(kfl_csysm_sld)
        call iexcha(kfl_prdef_sld)
        call iexcha(kfl_cohes_sld)
        do jr=1,3
           call iexcha(kfl_cshel_sld(jr))
        end do
        call iexcha(kfl_damag_sld)
        call iexcha(kfl_vofor_sld)
        call iexcha(kfl_plane_sld)
        call iexcha(kfl_restr_sld)
        call iexcha(kfl_indis_sld(1))
        call iexcha(kfl_indis_sld(2))
        call iexcha(kfl_invel_sld)
        call iexcha(kfl_plast_sld)
        do jr=1,3
           call iexcha(kfl_moduf_sld(jr))
        end do
        call iexcha(nmate_sld)
        do jr=1,3
           call iexcha(oripa_sld(jr))
        end do

        do jr=1,3
           call rexcha(gfibe_sld(jr))
           call rexcha(gravi_sld(jr))
           call rexcha(invel_sld(jr))
        end do
        call rexcha(grnor_sld)
        call rexcha(rbmas_sld)
        call rexcha(thick_sld)
        call rexcha(thiso_sld(1)) ! Displacement thershold
        call rexcha(thiso_sld(2)) ! Vonmises threshold
        do jr=1,9
           call rexcha(csysm_sld(jr))
        end do
        !
        ! Exchange of sld_reanut variables
        !
        call iexcha(kfl_xfeme_sld)
        call iexcha(kfl_xfcra_sld)
        call iexcha(kfl_stabi_sld)
        call iexcha(kfl_resid_sld)
        call iexcha(kfl_timet_sld)
        call iexcha(kfl_ninex_sld)
        call iexcha(kfl_tisch_sld)
        call iexcha(kfl_serei_sld)
        call iexcha(kfl_limit_sld)
        call iexcha(kfl_volca_sld)
        call iexcha(kfl_prest_sld)
        call iexcha(kfl_safet_table_sld)
        call iexcha(    miinn_sld)
        call iexcha(    minex_sld)
        call iexcha(kfl_plane_sld)
        call iexcha(kfl_gdepo)

        call rexcha(tifac_sld(1))
        call rexcha(tifac_sld(2))
        call rexcha(tifac_sld(3))
        call rexcha(tifac_sld(4))
        call rexcha(tifac_sld(5))

        call rexcha(cotol_sld)
        call rexcha(safet_sld)
        call rexcha(epsex_sld)
        call rexcha(safet_pseud_sld)
        call rexcha(factor_penal_sld)
        call rexcha(dafac_sld)
        call rexcha(safex_sld    )              ! Time function parameter for safety factor
        call rexcha(safma_sld    )              ! Maximum safety factor
        call rexcha(sstol_sld)
        call rexcha(masss_sld)

        call iexcha(mcavi_sld)
        call iexcha(nisaf_sld    )              ! Initial time step for variable cfl
        do ji=1,4
           call iexcha(kcavi_sld(  ji))
           call iexcha(iocav_sld(  ji))
           call rexcha(ocavi_sld(1,ji))
           call rexcha(ocavi_sld(2,ji))
           call rexcha(ocavi_sld(3,ji))
        end do
        !
        ! Exchange data read in sld_reaous
        !
        call iexcha(kfl_csysp_sld)
        call iexcha(kfl_exacs_sld)
        call iexcha(kfl_foten_sld)
        call iexcha(kfl_rotei_sld)

        call iexcha(kfl_rsfor_sld)
        call iexcha(kfl_pseud_sld)
        call iexcha(kfl_penal_sld)
        do ji=1,10
           call iexcha(idilimi_sld(ji))
           call rexcha(reslimi_sld(ji))
           call rexcha(safet_table_sld(1,ji))   ! Safety factor table: factor
           call rexcha(safet_table_sld(2,ji))   ! Safety factor table: residual
        end do

        call rexcha(rorig_sld(1))
        call rexcha(rorig_sld(2))
        call rexcha(rorig_sld(3))
        do ji = 1,10
           do jr=1,9
              call rexcha(csysp_sld(jr,ji))
           end do
        end do
        call posdef(1_ip,idumy)
        !
        ! Solvers
        !
        call soldef(1_ip)
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','sld_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','sld_sendat',parre)
           if(ISLAVE) call Parall(2_ip)
        end if
     end do

     if(IMASTER) call Parall(2_ip)

     call memchk(2_ip,istat,mem_servi(1:2,ID_PARALL),'parin','sld_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','sld_sendat',0_ip)
     call memchk(2_ip,istat,mem_servi(1:2,ID_PARALL),'parre','sld_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','sld_sendat',0_ip)

     !----------------------------------------------------------------
     !
     ! Exchange material properties read in sld_reaphy
     !
     !----------------------------------------------------------------

     !
     ! Allocatable related to materials
     !
     call Parall(30_ip)
     if( ISLAVE ) call sld_memphy(1_ip)

     do parii = 1,2
        npari = 0
        nparr = 0
        nparc = 0

        do ii = 1,nmate_sld
           call iexcha(modfi_sld(ii))
           call iexcha(modor_sld(1,ii))
           call iexcha(modor_sld(2,ii))
           call iexcha(kfl_dampi_sld(ii))
           call iexcha(kfl_coupt_sld(ii))
           call iexcha(lawco_sld(ii))
           call iexcha(lawst_sld(ii))
           call iexcha(lawho_sld(ii))
           call iexcha(lawch_sld(ii))
           call iexcha(lawpl_sld(ii))
           call iexcha(lawta_sld(ii))
           call rexcha(cocof_sld(ii))
           call rexcha(preti_sld(ii))
           call rexcha(timec_sld(ii))
           call rexcha(hillc_sld(ii))
           call rexcha(cal50_sld(ii))
           call rexcha(ortk1_sld(ii))
           call rexcha(ortk2_sld(ii))
           call rexcha(trans_sld(1,ii))
           call rexcha(trans_sld(2,ii))
           call iexcha(kfl_eccty(ii))  ! kernel variable, but read by solidz module
        end do

        do ii = 1,nmate_sld
           do ji = 1,ncoef_sld
              call rexcha(parsp_sld(ji,ii))
              call rexcha(densi_sld(ji,ii))
              call rexcha(parco_sld(ji,ii))
              call rexcha(parcc_sld(ji,ii))
              call rexcha(parch_sld(ji,ii))
              call rexcha(parcf_sld(ji,ii))
              call rexcha(velas_sld(ji,ii))
           end do
           call rexcha(dampi_sld(1,ii))
           call rexcha(dampi_sld(2,ii))
        end do
        !
        ! Undamaged stiffness tensor                     (AQU)
        !
        do ii = 1,nmate_sld
           do ji = 1,ndime*2
              do ki = 1,ndime*2
                 call rexcha(stiff0_sld(ji,ki,ii))
              end do
           end do
        end do

        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,5),'parin','sld_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,5),'parre','sld_sendat',parre)
           if( ISLAVE ) call Parall(two)
        end if
     end do

     if( IMASTER ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','sld_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','sld_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','sld_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','sld_sendat',0_ip)

     !-------------------------------------------------------------------
     !
     ! Dimensions read in sld_reabcs
     !
     !-------------------------------------------------------------------

     call Parall(27_ip)

     strin = 'SLD_REABCS'
     strre = 'SLD_REABCS'

     do parii = 1,2
        npari = 0
        nparr = 0
        nparc = 0

        call iexcha(kfl_conbc_sld)
        call iexcha(kfl_local_sld)
        call iexcha(kfl_csysl_sld)
        call iexcha(kfl_cycle_sld)
        call iexcha(nfunc_sld)
        call iexcha(ncrak_sld)
        call iexcha(kfl_bodyf_sld)
        call iexcha(kfl_windk_sld)
        call iexcha(kfl_conta_sld)
        call iexcha(contactbou_sld)
        call iexcha(coupling_contact_its)

        call rexcha(neumann_relax)
        call rexcha(coupling_contact_tol)
        do ki =1,3
           call rexcha(vect_proje_sld(ki))
        end do
        do ki =1,9
           call rexcha(csysl_sld(ki))
        end do

        do ki =1,4
           call rexcha(tzero_sld(ki))
           call rexcha(tpstr_sld(ki))
           call rexcha(pzero_sld(ki))
           call rexcha(pstr0_sld(ki))
           call rexcha(cpres_sld(ki))
           call rexcha(rpres_sld(ki))
           call rexcha(part0_sld(ki))
           call rexcha(ppost_sld(ki))
           call rexcha(gfill_sld(ki))
           call rexcha(prest_sld(ki))
        end do

        do ki = 1,10
           call iexcha(mtloa_sld(ki))
        end do
        do ki = 1,20
           call rexcha(rtico_sld(1,ki))
           call rexcha(rtico_sld(2,ki))
           do ji = 1,20
              call rexcha(fubcs_sld(ki,ji))
              call iexcha(kfl_funty_sld(ki,ji))
           end do
        end do

        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','sld_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','sld_sendat',parre)
           if( ISLAVE ) call Parall(2_ip)
        end if
     end do

     if( IMASTER ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','sld_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','sld_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','sld_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','sld_sendat',0_ip)

     !-------------------------------------------------------------------
     !
     ! Arrays read in sld_reabcs
     !
     !-------------------------------------------------------------------

     if( ISLAVE .and. ncrak_sld > 0 ) call sld_membcs(31_ip) ! Cracks
     call spnbcs(tncod_sld)
     !call spnbcs(tgcod_sld)
     call spbbcs(tbcod_sld)

     if( ISLAVE ) then
        call sld_membcs(3_ip)
        do ifunc=1,10
           if (mtloa_sld(ifunc) == 0) mtloa_sld(ifunc) = 1            ! done to allocate a default memory space
           call sld_membcs(10_ip + ifunc)  ! allocate the prescription time function vector for ifunc
        end do
     end if

     if( ndime == 2 ) then
        idumy = 2
     else
        idumy = 4
     end if

     do parii=1,2
        npari=0
        nparr=0
        nparc=0

        do ifunc=1,10
           do ki=1,mtloa_sld(ifunc)
              do ji= 1,ndime+1
                 call rexcha(tload_sld(ifunc)%a(ji,ki))
              end do
           end do
        end do
        do ki = 1,ncrak_sld
           do ii = 1,idumy
              do ji = 1,ndime
                 call rexcha(crkco_sld(ji,ii,ki))
              end do
           end do
        end do

        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','sld_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','sld_sendat',parre)
           if( ISLAVE ) call Parall(2_ip)
        end if
     end do

     if( IMASTER ) call Parall(2_ip)
     if( IMASTER .and. ncrak_sld > 0 ) call sld_membcs(-31_ip) ! Cracks

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','sld_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','sld_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','sld_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','sld_sendat',0_ip)

  end select

  npari = 0
  nparr = 0
  nparc = 0

end subroutine sld_sendat
