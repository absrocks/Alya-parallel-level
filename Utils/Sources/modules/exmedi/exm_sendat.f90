subroutine exm_sendat(order)
  !-----------------------------------------------------------------------
  !****f* exmedi/exm_sendat
  ! NAME
  !    exm_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    exm_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_inpout
  use      mod_memchk
  use      mod_opebcs
  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,kr,idumy
  integer(ip)             :: imate,iroot,idime
  integer(ip)             :: ibcas,ixchn,kfl_ptask_old 
  integer(4)              :: istat

  ibcas= 2_ip
  ixchn= 300_ip


  select case (order)

  case(1)
     !
     ! Exchange data read in exm_reaphy, exm_reanut and exm_reaous
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of exm_reaphy variables 
        !        
        call iexcha(kfl_timei_exm)
        call iexcha(kfl_gemod_exm)
        call iexcha(kfl_fento_exm)
        call iexcha(kfl_appli_exm)
        call iexcha(kfl_appty_exm)
        call iexcha(kfl_appva_exm)
        call iexcha(kfl_gcoup_exm)
        call iexcha(kfl_comdi_exm)
        call iexcha(kfl_stead_exm)
        call iexcha(kfl_cemod_exm)
        call iexcha(kfl_ptrig_exm)
        call iexcha(kfl_drugs_exm)
        call iexcha(kfl_prdef_exm)
        call iexcha(kfl_stree_exm)
        call iexcha(kfl_heter_exm)
        call iexcha(kfl_nodif_exm)
        call iexcha(kfl_atbhe_exm)
        call iexcha(kfl_paced_exm)
        call iexcha(onecl_exm(1))               ! Number of beats to steady state
        call iexcha(onecl_exm(2))               ! Heart rate to set initial conditions
        call iexcha(ndofn_exm)
        call iexcha(ndof2_exm)
        call iexcha(nmate_exm)
        call iexcha(nevat_exm)
        call iexcha(ncomp_exm)

        call iexcha(nrootecg_exm)
        do iroot = 1,256_ip                     ! All possible lead roots
           do idime = 1,ndime
              call rexcha(coordecg_exm(idime,iroot))
           end do
        end do
        call iexcha(kfl_psecg_exm)
        call iexcha(kcopeecg_exm)
        call iexcha(nrootecg_exm)
        call iexcha(nvint_exm)
        
        do imate=1,nmate
           call iexcha(kfl_hetermate_exm(imate))
           call iexcha(kfl_fract_diffusion_exm(imate))
           call iexcha(kfl_inaga_exm(imate))
           call iexcha(kfl_hfmod_exm(imate))
           call iexcha(kfl_hfmodmate_exm(imate))
           call rexcha(xmccmmate_exm(imate))
           call iexcha(moneclmate_exm(1,imate))
           call iexcha(moneclmate_exm(2,imate))

           call iexcha(kfl_drugsmate_exm(imate))

           do jr=1,12
              call rexcha(ttparmate_exm(1,jr,imate))
              call rexcha(ttparmate_exm(2,jr,imate))
              call rexcha(ttparmate_exm(3,jr,imate))
              call rexcha(drugdmate_exm(  jr,imate))
           end do

        end do

        call iexcha(ngrou_exm)
        call iexcha(nauxi_exm)
        call iexcha(nconc_exm)
        call iexcha(nicel_exm)
        call iexcha(nstim_exm)
        call iexcha(nstis_exm)
        call iexcha(modfi_exm)
        call iexcha(modac_exm)
        call iexcha(modce_exm)
        call iexcha(modab_exm)
        call iexcha(modst_exm)

        call rexcha(dtinv_exm)
        call rexcha(cleng_exm)                  ! Reference characteristic length
        call rexcha(scond_exm)
        call rexcha(react_exm)
        call rexcha(xthri_exm)
        call rexcha(xthrs_exm)

        call rexcha(thiso_exm(1))               ! Isochrones trigger
        call rexcha(thiso_exm(2))               ! Isochrones threshold
        call rexcha(thiso_exm(3))               ! Isochrones auto
        do jr=1,3
           do kr=1,12
              call rexcha(ttpar_exm(jr,kr))                  ! Current Remodelling of cell models
           end do
        end do
        do jr=1,12
           call rexcha(drugd_exm(jr))                  !  Drug definitions
        end do

        do jr= 1,15000

           call rexcha(aptim(jr))
           call rexcha(apval_exm(jr))
           call rexcha(aplap_exm(jr))
           call rexcha(aprea_exm(jr))
           do kr=1,3
              call rexcha(apcen_exm(kr,jr))
           end do
        end do

        do jr=1,2
           do kr=1,2
              call rexcha(ceglo_exm(jr,kr))
              do imate=1,nmate
                 call rexcha(gdiff_exm(jr,kr,imate))
              end do
           end do
        end do

        call rexcha(gfibe_exm(1))
        call rexcha(gfibe_exm(2))
        call rexcha(gfibe_exm(3))

        !
        ! Streeter fiber model
        !
        call iexcha(nstrb_exm   )
        call rexcha(strbo_exm(1))
        call rexcha(strbo_exm(2))
        call rexcha(strbo_exm(3))
        call rexcha(fiaxe_exm(1))
        call rexcha(fiaxe_exm(2))
        call rexcha(fiaxe_exm(3))
        call rexcha(stran_endo_exm)
        call rexcha(stran_epi_exm)


        do imate= 1,nmate
           call iexcha(kfl_cellmod(imate))
        end do

        do imate=1,nmate
           do jr=1,20
              call rexcha(xmopa_exm(jr,imate))
           end do
           call rexcha(voini_exm(imate))
           call iexcha(kfl_voini_exm(imate))          
           call iexcha(kfl_fract_diffusion_exm(imate))
           call rexcha(fract_diff_exm(imate))
        end do


        call rexcha(xmccm_exm)
        call rexcha(aploo_exm(1))
        call rexcha(aploo_exm(2))
        call rexcha(aploo_exm(3))     
        call rexcha(poref_fhn_exm(1))
        call rexcha(poref_fhn_exm(2))

        !
        ! Exchange of exm_reanut variables 
        !

        call iexcha(kfl_genal_exm)              ! General alg. type
        call iexcha(kfl_goite_exm)              ! Keep iterating


        call iexcha(kfl_shock_exm)              ! 
        call iexcha(kfl_comat_exm)              ! 
        call iexcha(kfl_weigh_exm)              ! 
        call iexcha(kfl_timet_exm)              ! 
        call iexcha(kfl_tiacc_exm)              ! 
        call iexcha(kfl_ticel_exm)              ! 
        call iexcha(kfl_tisch_exm)              ! 
        call iexcha(kfl_normc_exm)              ! 
        call iexcha(kfl_adres_exm)              ! 
        call iexcha(kfl_algso_exm)              ! 
        call iexcha(kfl_repro_exm)              ! 
        call iexcha(kfl_nolim_exm)              ! 
        call iexcha(kfl_nolum_exm)              ! 
        call iexcha(miinn_exm)              ! 
        call iexcha(msste_exm)              ! 
        call iexcha(mnoli_exm)              ! 
        call iexcha(msoit_exm)              ! 
        call iexcha(nkryd_exm)              ! 
        call iexcha(itera_exm)              ! 
        call iexcha(nunkn_exm)              ! 


        call rexcha(dtcri_exm    )              ! 
        call rexcha(shock_exm    )              ! 
        call rexcha(sstol_exm    )              ! 
        call rexcha(cotol_exm    )              ! 
        call rexcha(corat_exm    )              ! 
        call rexcha(dtext_exm    )              ! 
        call rexcha(safet_exm    )              ! 
        call rexcha(solco_exm    )              ! 
        call rexcha(weigh_exm    )              ! 
        call rexcha(tnoli_exm    )              ! 

        do ji=1,2
           call rexcha(err01_exm(ji))           ! L1 error u
           call rexcha(err02_exm(ji))           ! L2 error u
           call rexcha(err0i_exm(ji))           ! Linf error u
           call rexcha(err11_exm(ji))           ! L1 error grad(u)
           call rexcha(err12_exm(ji))           ! L2 error grad(u)
           call rexcha(err1i_exm(ji))           ! Linf error grad(u)
           call rexcha(cpu_exmed(ji))           ! CPU for the exm problem
        end do
        do ji=1,3
           call rexcha(staco_exm(ji))           ! 
        end do
        do ji=1,10
           call rexcha(resid_exm(ji))           ! Residual for outer iterations (Mom,Cont,Ene,Glo)
        end do

        !
        ! Exchange data read in exm_reaous
        !

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
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','exm_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','exm_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(ibcas)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(ibcas)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','exm_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','exm_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','exm_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','exm_sendat',0_ip)



     !     if(nmate_exm>=2) then
     !        
!!!!        if(kfl_paral>=1) call exm_memphy(1_ip)
     !
     !        call runend('EXM_SENDAT: NOT PREPARED YET')
     !
     !      party = 1  --> vector dimensioned nelem
     !      party = 2  --> vector dimensioned nboun
     !      party = 3  --> vector dimensioned npoin
     !      party = 4  --> vector dimensioned nbopo
     !
     !      pardi = 1, 2, 3 -->  number of columns of a vector (pari1, 2,3 )
     !      pard1 = ...     -->  size of the first column
     !      parki = 1, 2, 3 -->  integer, real or character
     !

     !     end if
     kfl_ptask = kfl_ptask_old

  case (2)

     !
     ! Exchange data read in exm_reabcs

     !
     !
     ! Allocate memory for the first pass
     !

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0

        call iexcha(kfl_exboc_exm)

        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','exm_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','exm_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(ibcas)
        end if

     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(ibcas)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','exm_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','exm_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','exm_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','exm_sendat',0_ip)

     if(kfl_paral/=0.or.kfl_ptask/=2) then

        if (kfl_exboc_exm == 1) then

           call spnbcs(tncod_exm)
           call spbbcs(tbcod_exm)

        end if

     end if

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine exm_sendat
