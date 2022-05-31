
subroutine nsa_sendat(order)
  !-----------------------------------------------------------------------
  !****f* nastal/nsa_sendat
  ! NAME
  !    nsa_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    nsa_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  use      def_inpout
  use      def_solver
  use      mod_memchk
  use      mod_opebcs

  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,jc,ki,kr,kc,idumy
  integer(ip)             :: ielem,ipoin
  integer(ip)             :: jnode,jpoin,jelem,inode
  integer(ip)             :: kfl_ptask_old 
  integer(4)              :: istat,ifunc

  select case (order)

  case(1)     
     !
     ! Exchange data read in nsa_reaphy, nsa_reanut and nsa_reaous
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
      

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !-------------------------------------------------------------------------
        ! Exchange of nsa_reaphy variables 
        !-------------------------------------------------------------------------

        call iexcha(kfl_stead_nsa)
        call iexcha(kfl_advec_nsa)
        call iexcha(kfl_visco_nsa)
        call iexcha(kfl_inico_nsa)
        call iexcha(kfl_infun_nsa)
        call iexcha(kfl_benme_nsa)
        call iexcha(kfl_botta_nsa)
        call iexcha(kfl_rearst_nsa)
        call iexcha(kfl_sptyp_nsa)
        call iexcha(kfl_hysta_nsa)
        call iexcha(kfl_inleb_nsa)
        call iexcha(kfl_dampf_nsa)
        call iexcha(kfl_mfrco_nsa)
        call iexcha(kfl_brunt_nsa)
        call iexcha(kfl_coupl_nsa)
        call iexcha(kfl_cotur_nsa)
        call iexcha(kfl_sponge_nsa)
        call iexcha(kfl_cylind_nsa)
        call iexcha(kfl_refpdt_nsa)
        
        !------------------------------------------------------------------------
        ! Additional INTEGER variables used only locally with METEO problems:
        !------------------------------------------------------------------------
        call iexcha(kfl_adiff_nsa)
        call iexcha(vtkrestart_nsa)
        call iexcha(istep_nsa)
        call iexcha(nvar_nsa)
        call iexcha(nelx_nsa)
        call iexcha(nely_nsa)
        call iexcha(nelz_nsa)
        call iexcha(nx_nsa)
        call iexcha(ny_nsa)
        call iexcha(nz_nsa)
        call iexcha(ncol_nsa)
        call iexcha(kfl_sptyp_nsa)

        call iexcha(kfl_uniformvelocity_nsa)
        call iexcha(kfl_sponge_nsa)
        call iexcha(kfl_thetac_nsa)
        call iexcha(kfl_tracerc_nsa)
        call iexcha(kfl_xc_nsa)
        call iexcha(kfl_yc_nsa)
        call iexcha(kfl_zc_nsa)
        call iexcha(kfl_rc_nsa) 
        call iexcha(kfl_xr_nsa)
        call iexcha(kfl_yr_nsa)
        call iexcha(kfl_zr_nsa)
        call iexcha(kfl_rc_nsa)
        call iexcha(kfl_spher_nsa)
        !------------------------------------------------------------------------
        ! END additional INTEGER variables used only locally with METEO problems.
        !------------------------------------------------------------------------

        call iexcha(nmate_nsa)
        call iexcha(ndofn_nsa)
        call iexcha(nkeep_nsa)
        call iexcha(ndof2_nsa)
        call iexcha(nevat_nsa)
        call iexcha(nevab_nsa)
        call iexcha(nflub_nsa)
        call iexcha(ivert_nsa)
        call iexcha(lawde_nsa)
        call iexcha(lawvi_nsa)
        call iexcha(lawpo_nsa)
        call iexcha(lawst_nsa)
        call iexcha(mfrse_nsa)
        call iexcha(kfl_inibu_nsa)              ! Built-in initial fields type
        call iexcha(kfl_inifi_nsa(1))           ! Initial fields 
        call iexcha(kfl_inifi_nsa(2))           ! Initial fields 
        call iexcha(kfl_inifi_nsa(3))           ! Initial fields 
        call iexcha(kfl_inkee_nsa(2))           ! Hydrostatic reference, density
        call iexcha(kfl_inkee_nsa(3))           ! Hydrostatic reference, temperature
        call iexcha(kfl_inkee_nsa(4))           ! Hydrostatic reference, pressure

        call iexcha(kfl_iposi_nsa   )           ! Positional initial fields 
        ! Positional initial fields definition
        do jr=1,8
           call iexcha(iposi_nsa(jr)%kflag)
           call iexcha(iposi_nsa(jr)%geometry)
           call rexcha(iposi_nsa(jr)%center(1))
           call rexcha(iposi_nsa(jr)%center(2))
           call rexcha(iposi_nsa(jr)%center(3))
           call rexcha(iposi_nsa(jr)%radius)        
           call rexcha(iposi_nsa(jr)%value)        
        end do

        do jr=1,10
           call iexcha(nfiel_nsa(jr))           ! Fields assignement
        end do
        
        call rexcha(dtinv_nsa)
        call rexcha(dtlim_nsa)
        call rexcha(adgam_nsa)                  ! Adiabatic exponent
        call rexcha(cpcoe_nsa)                  ! Cp coefficient
        call rexcha(cvcoe_nsa)                  ! Cv coefficient
        call rexcha(rgasc_nsa)                  ! R gas constant
        call rexcha(rgacv_nsa)                  ! R gas constant over Cv
        call rexcha(fcons_nsa)                  ! Convection term
        call rexcha(fvins_nsa)                  ! Viscous term
        call rexcha(entro_nsa)                  ! Reference inflow entropy
        call rexcha(afact_nsa)                  ! Co factor for potential temperature
        call rexcha(tstag_nsa)                  ! Stagnation temperature
        call rexcha(trefa_nsa)                  ! Temperature recovery factor
        call rexcha(cppra_nsa)                  ! Cp / Prandtl
        call rexcha(cpprt_nsa)                  ! Cp / Turbulent Prandtl
        call rexcha(pbaro_nsa)                  ! Reference pressure at ground level for barometric law
        call rexcha(sofac_nsa)                  ! Sound speed factor
        call rexcha(turbu_nsa)                  ! Eddy viscosity coefficient
        call rexcha(densi_nsa)                  ! Reference density (rho)
        call rexcha(mfrgr_nsa)                  ! Mass flow rate growth parameter
        call rexcha(mfrub_nsa)                  ! Bulk velocity Ub target
        call rexcha(ubpre_nsa)                  ! Bulk velocity Ub target previous time-step
        call rexcha(tempe_nsa)                  ! Reference temperature (rho)
        call rexcha(mowei_nsa)                  ! Molecular weight of the mixture
        call rexcha(runiv_nsa)                  ! Molecular weight of the mixture

        do jr=1,ndofn_nsa
           do kr=1,ndofn_nsa
              call rexcha(gravm_nsa(kr,jr))     ! Gravity matrix
           end do
        end do
   
        call rexcha(mixrv_nsa)                  ! Vapor mixing ratio
         
        call rexcha(speed_nsa)                  ! Reference velocity module (u)
        call rexcha(press_dynamic_nsa)          ! Reference pressure (p)
        call rexcha(press_nsa)                  ! Reference pressure (p)
        call rexcha(visco_nsa)                  ! Reference dynamic viscosity (mu)
        call rexcha(thdif_nsa)                  ! Reference thermal diffusion (k)
        call rexcha(prand_nsa)                  ! Prandtl number (SUPERSEDES k)
        call rexcha(pratu_nsa)                  ! Turbulent Prandtl number
        call rexcha(rreyn_nsa)                  ! Reference Reynolds number (SUPERSEDES mu)
        call rexcha(rmach_nsa)                  ! Reference Mach number (SUPERSEDES Cp)
        call rexcha(attxy_nsa)                  ! Attack angle plane X-Y
        call rexcha(attxz_nsa)                  ! Attack angle plane X-Z
        call rexcha(attyz_nsa)                  ! Attack angle plane Y-Z
        call rexcha(axyin_nsa)                  ! Attack angle plane X-Y (for inlet boundary conds) 
        call rexcha(axzin_nsa)                  ! Attack angle plane X-Z (for inlet boundary conds) 
        call rexcha(ayzin_nsa)                  ! Attack angle plane Y-Z (for inlet boundary conds) 
        call rexcha(spein_nsa)                  ! Velocity module (for inlet boundary conds)
        call rexcha(cleng_nsa)                  ! Reference characteristic length
        call rexcha(grnor_nsa)
        call rexcha(brure_nsa)                  ! Brunt frequency reference value
        !More meteo:
        call rexcha(xmin_nsa)
        call rexcha(ymin_nsa)
        call rexcha(zmin_nsa)
        call rexcha(xmax_nsa)
        call rexcha(ymax_nsa)
        call rexcha(zmax_nsa)
        call rexcha(uvelo_nsa)
        call rexcha(vvelo_nsa)
        call rexcha(wvelo_nsa)
        call rexcha(xradi_nsa)
        call rexcha(yradi_nsa)
        call rexcha(zradi_nsa)
        call rexcha(xc_nsa)
        call rexcha(yc_nsa)
        call rexcha(zc_nsa)
        call rexcha(dxs_nsa)
        call rexcha(dys_nsa)
        call rexcha(dzs_nsa)
        call rexcha(ampx_nsa)
        call rexcha(ampy_nsa)
        call rexcha(ampz_nsa)
        call rexcha(thetac_nsa)
        call rexcha(tracerc_nsa)

        !built in models: laxliu
        do jr=1,4
           call rexcha(llrho_nsa(jr))
           call rexcha(llpre_nsa(jr))
           call rexcha(llvel_nsa(1,jr))
           call rexcha(llvel_nsa(2,jr))
        end do
        call rexcha(llcen_nsa(1))
        call rexcha(llcen_nsa(2))
        call rexcha(llcen_nsa(3))

        do jr=1,ncoef_nsa
           call rexcha(vispa_nsa(jr))           ! Parameters for the law of viscosity
        end do
        do jr=1,ncoef_nsa
           call rexcha(poros_nsa(jr))           ! Porosity (sig)
        end do
        do jr=1,3
           call rexcha(veloc_nsa(jr))           ! Reference velocity (rho)
           call rexcha(gravi_nsa(jr))
           call rexcha(vranr_nsa(jr,1))
           call rexcha(vranr_nsa(jr,2))
           call rexcha(vranr_nsa(jr,3))
           call rexcha(frayl_nsa(jr,1))
           call rexcha(frayl_nsa(jr,2))
           call rexcha(frayl_nsa(jr,3))
           call rexcha(frayl_nsa(jr,4))
        end do
        do jr=1,5
           call rexcha(teano_nsa(jr))
           call rexcha(xbubb_nsa(1,jr))
           call rexcha(xbubb_nsa(2,jr))
           call rexcha(xbubb_nsa(3,jr))
           call rexcha(xbubb_nsa(4,jr))
           call rexcha(xbubb_nsa(5,jr))
           call rexcha(xbubb_nsa(6,jr))
        end do
        do jr=1,5
           call rexcha(rgcou_nsa(jr))
           call rexcha(cpcou_nsa(jr))
           call rexcha(cvcou_nsa(jr))
        end do

        !
        ! Exchange of nsa_reanut variables 
        !
        call iexcha(kfl_mod_elmop_nsa)
        call iexcha(kfl_timei_nsa)
        call iexcha(kfl_algor_nsa)              ! Algorithm tipe (monolithic, fractional...)
        call iexcha(kfl_timet_nsa)              ! Explicit / Implicit
        call iexcha(kfl_matri_nsa)              ! Explicit with matrix
        call iexcha(kfl_diagi_nsa)              ! Diagonal implicit terms
        call iexcha(kfl_lotim_nsa)              ! Diagonal implicit terms
        call iexcha(kfl_foreg_nsa)              ! Forced regime: compressible or incompressible
        call iexcha(kfl_relat_nsa)              ! Relativistic flow
        call iexcha(kfl_isent_nsa)              ! Isentropic flow
        call iexcha(kfl_ncons_nsa)              ! Non-conservative set
        call iexcha(kfl_unkse_nsa)              ! Conservative or primitive unknowns set
        call iexcha(kfl_rayle_nsa)              ! Rayleigh dumping
        call iexcha(kfl_lopre_nsa)              ! Local preconditioning or not
        call iexcha(kfl_zero_initial_velocity_nsa) ! Identifies those problems with zero initial velocity
        call iexcha(kfl_pseud_nsa)              ! Pseudo time step is used or not
        call iexcha(kfl_locti_nsa)              ! Local or global time step
        call iexcha(kfl_goite_nsa)              ! Keep iterating
        call iexcha(kfl_stabi_nsa)              ! Stabilization method
        call iexcha(kfl_galer_nsa)              ! Galerkin terms form
        call iexcha(kfl_repro_nsa)              ! Stabilization based on residual projection
        call iexcha(kfl_hconv_nsa)              ! Convective characteristic length
        call iexcha(kfl_stafl_nsa)              ! 
        call iexcha(kfl_penal_nsa)              ! Penalization
        call iexcha(kfl_turbu_nsa)              ! Turbulent flow
        call iexcha(kfl_weigh_nsa)              ! Weighting of du/dt
        call iexcha(kfl_tiacc_nsa)              ! Temporal accuracy
        call iexcha(kfl_normc_nsa)              ! Norm of convergence
        call iexcha(kfl_taudi_nsa)              ! Diagonal tau
        call iexcha(kfl_track_nsa)              ! Subscales tracking
        call iexcha(kfl_higha_nsa)              ! High aspect ratio elements are present
        call iexcha(kfl_resmo_nsa)              ! Residual smoothing
        call iexcha(kfl_skews_nsa)              ! 
        call iexcha(kfl_linea_nsa)              ! Linearization (RHS=0, Picard=1, Newton=2)
        call iexcha(kfl_delun_nsa)              ! Delta form
        call iexcha(kfl_algso_nsa)              ! Type of algebraic solver
        call iexcha(kfl_tisch_nsa)              ! Time integration scheme
        call iexcha(kfl_tiext_nsa)              ! Externally fixed time step flag
        call iexcha(kfl_cfllo_nsa)              ! Local CFL
        call iexcha(kfl_dtadj_nsa)              ! Time increment type (global or local) adjustment
        call iexcha(kfl_fasts_nsa)              ! Fast stationary (number of neighbors layers)
        call iexcha(kfl_nacdr_nsa)              ! Nastal cdr
        call iexcha(kfl_reate_nsa)              ! Reate in tau and dt
        call iexcha(miinn_nsa)                  ! Max # of iterations
        call iexcha(minew_nsa)                  ! Initial Newton iteration
        call iexcha(mcfll_nsa)                  ! Max # of CFL values for the CFL_LOCAL option
        call iexcha(miinn_pseud_nsa)                  ! Max # of pseudo-time iterations
        call iexcha(msoit_nsa)                  ! Max # of solver iterations
        call iexcha(npica_nsa)                  ! Number of Picard iteration (Newton's lin.)
        call iexcha(nkryd_nsa)                  ! Krylov dimension
        call iexcha(neule_nsa)                  ! # of Euler time steps
        call iexcha(nunkn_nsa)                  ! # of unknonws = ndofn_nsa*npoin  
        do ji=1,10
           call iexcha(kfl_dttyp_nsa(ji))       ! Time increment type (global or local) per equation
        end do
        do ji=1,4
           call iexcha(kfl_ximpl_nsa(ji))       ! Explicit terms in the implicit form
        end do
        do ji=1,10
           call iexcha(kfl_taufa_nsa(ji,1))     ! Time factors
           call iexcha(kfl_taufa_nsa(ji,2))     ! Tau factors
        end do
        do ji=1,5
           call iexcha(kfl_shock_nsa(ji))       ! Shock capturing type 
        end do
        do ji=1,3
           call iexcha(kranr_nsa(ji,1))           ! 
           call iexcha(kranr_nsa(ji,2))           ! 
           call iexcha(kranr_nsa(ji,3))           ! 
           call iexcha(kranr_nsa(ji,4))           ! 
        end do
        call iexcha(ncomp_nsa)
        call iexcha(nromp_nsa)
        call iexcha(nfrap_nsa)
        call iexcha(ntomp_nsa)
        do ji=1,500
           call rexcha(fcfll_nsa(ji,1))
           call rexcha(fcfll_nsa(ji,2))
        end do

    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
    call iexcha( euler_nsa )
    call iexcha( lodi_nsa  )
    call rexcha( prefe_lodi_nsa )
    call rexcha( sigma_lodi_nsa )
    !-----------------------------------------------------------------------! 

    call rexcha(shock_nsa    )              ! Shock capturing parameter
    call rexcha(shtol_nsa    )              ! Shock capturing residual tolerance
    call rexcha(sstol_nsa    )              ! Steady state tolerance
    call rexcha(cotol_nsa    )              ! Convergence tolerance
    call rexcha(corat_nsa    )              ! Convergence tolerance ratio
    call rexcha(penal_nsa    )              ! Penalization factor
    call rexcha(dtext_nsa    )              ! Externally fixed time step 
    call rexcha(safet_nsa    )              ! Safety factor for time step
    call rexcha(safex_nsa    )              ! Time function parameter for safety factor
    call rexcha(safma_nsa    )              ! Maximum safety factor
    call rexcha(safet_pseud_nsa)            ! Safety factor for pseudo time step
    call rexcha(safrk_nsa    )              ! Delta time shrinking factor for RK schemes
    call rexcha(solco_nsa    )              ! Solver tolerance
    call rexcha(weigh_nsa    )              ! Weight of dU/dt in the residual
    call rexcha(xfree_nsa    )              ! X Coordinate of the plane where to free 
    do ji=1,10
       call rexcha(dtcri_nsa(ji))           ! Critical time step
       call rexcha(dtmax_nsa(ji))           ! Maximum of the local time steps
       call rexcha(theta_nsa(ji))           ! "Explicitness" weighting factors
       call rexcha(safet_table_nsa(1,ji))   ! Safety factor table: factor
       call rexcha(safet_table_nsa(2,ji))   ! Safety factor table: residual
    end do

    call iexcha(nisaf_nsa    )              ! Initial time step for variable cfl 
    call iexcha(kfl_safet_table_nsa )       ! Safety factor table flag

!BEGIN NEW PART
        !
        ! Exchange data read in nsa_reabcs
        !
         

        call iexcha(kfl_modfi_nsa)              ! Modify kfl_fixno_nsa
        call iexcha(kfl_adres_nsa)              ! Adaptive residual tolerance for subiterations
        call iexcha(kfl_confi_nsa)              ! Confined flow
        call iexcha(kfl_local_nsa)              ! Local system of reference
        call iexcha(kfl_inlet_nsa)              ! Fix mass flow rate b.c.
        call iexcha(kfl_conbc_nsa)              ! Constant b.c.
        call iexcha(kfl_dumbo_nsa)              ! Dump derived boundary conditions
        call iexcha(kfl_spong_nsa   )           ! Rayleigh sponge flag
        call iexcha(kfl_skewb_nsa   )           ! Skew symetrics in bvess
        call iexcha(kfl_bofty_nsa)              ! Boundary conditions file type (nstinc or nastal)
        call iexcha(kfl_bofie_nsa)              ! Boundary conditions from fields, when present
        call iexcha(kfl_nopen_nsa)              ! No-penetration
        call iexcha(kfl_tredg_nsa)              ! Trailing edge condition
        call iexcha(nodpr_nsa)                  ! Node on which pressure is prescribed
        call iexcha(nfunc_nsa)                  ! Number of time-dependent b.c. functions

        !
        ! KFL_FUNTYP
        !
        do ki=1,20
           do ji=1,20
              call iexcha(kfl_funty_nsa(ki,ji))
           end do
        end do
        !
        ! FUNPA
        !
        do ki=1,20
           do ji=1,20
              call rexcha(fubcs_nsa(ki,ji))      
           end do
        end do
        !
        ! RTICO
        !
        do ji=1,20
           call rexcha(rtico_nsa(ji,1))      
           call rexcha(rtico_nsa(ji,2))      
        end do

        call rexcha(delta_nsa)                   ! Distance to the wall
        call rexcha(angle_tredg_nsa)             ! Angle threshold for the trailing edge detection (DEGREES)
        call rexcha(inlet_nsa)                   ! Target mass flow rate at the boundary
        
        !Time dependent boundaries
        do ki=1,10
           call iexcha(mtloa_nsa(ki))
        end do

!END NEW PART
        !
        ! Exchange data read in nsa_reaous
        !

        call posdef(1_ip,idumy)
       
        call iexcha(kfl_exacs_nsa)              ! Exact solution for the NS eq.
        do ji=1,2
           call iexcha(kfl_chkpo_nsa(ji))       ! Checkpoint-restart file flags: input(1) and output(2)
        end do
        call iexcha(kfl_bodyf_nsa)              ! Body force
        do ji=1,10
           do ki=1,2
              call rexcha(bodyr_nsa(ji,ki,1))
              call rexcha(bodyr_nsa(ji,ki,2))
              call rexcha(bodyr_nsa(ji,ki,3))
           end do
        end do
        !
        ! Solvers
        !
        call soldef(1_ip)
        !
        ! Sets
        !
        nparc=nparc+150
        if(parii==2.and.kfl_paral==0) parch(1:150)     = fil_chkpo_nsa(1)
        if(parii==2.and.kfl_paral>=1) fil_chkpo_nsa(1) = parch(1:150)       
        nparc=nparc+150
        if(parii==2.and.kfl_paral==0) parch(1:150)     = fil_chkpo_nsa(2)
        if(parii==2.and.kfl_paral>=1) fil_chkpo_nsa(2) = parch(1:150)  
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','nsa_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','nsa_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call par_broadc()
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'parin','nsa_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','nsa_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','nsa_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','nsa_sendat',0_ip)
     
     !
     ! Allocatable arrays: LMATE_NSA and LMATN_NSA
     !

      

     if(nmate_nsa>=2) then

        if(kfl_paral>=1) call nsa_memphy()

        call vocabu(111_ip,0_ip,0_ip)
        pari1 => lmate_nsa           
        call par_mygather()

        call vocabu(311_ip,0_ip,0_ip)
        pari1 => lmatn_nsa           
        call par_mygather()       

     end if
     kfl_ptask = kfl_ptask_old

     !-------------------------------------------------------------------
     !
     ! Arrays read in sld_reabcs
     !
     !-------------------------------------------------------------------
     call spnbcs(tncod_nsa)
     call spbbcs(tbcod_nsa)

     if( ISLAVE ) then
        call nsa_membcs(4_ip) 
        do ifunc=1,10
           if (mtloa_nsa(ifunc) == 0) mtloa_nsa(ifunc) = 1            ! done to allocate a default memory space
           do ji=1,mtloa_nsa(ifunc)
              call nsa_membcs(10_ip + ifunc)  ! allocate the prescription time function vector for ifunc
           end do
        end do
     end if

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        
        do ifunc=1,10
           do ki=1,mtloa_nsa(ifunc)
              do ji= 1,ndofn_nsa+1
                 call rexcha(tload_nsa(ifunc)%a(ji,ki))      
              end do
           end do
        end do
        
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','nsa_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','nsa_sendat',parre)
           if( ISLAVE ) call par_broadc()
        end if
     end do
     
     if( IMASTER ) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'parin','nsa_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','nsa_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','nsa_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','nsa_sendat',0_ip)


  case (2)   !! DMM this has to be removed, not used anymore!

     strin = 'NSA_REABCS'
     strre = 'NSA_REABCS'

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Allocate memory for the first pass
        !

        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','nsa_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','nsa_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call par_broadc()
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'parin','nsa_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','nsa_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','nsa_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','nsa_sendat',0_ip)

     !
     ! Open memory for boundary conditions
     !
     if(kfl_paral>=1) then
        call nsa_membcs(one)
        call nsa_membcs(two)
     end if

     if(kfl_paral/=0.or.kfl_ptask/=2) then
        !
        ! KFL_FIXNO_NSA
        !

        call vocabu(NPOIN_INTE_2DIM,ndofn_nsa,0_ip)
        pari2 => kfl_fixno_nsa
        call par_mygather()
        !
        ! KFL_FIXBO_NSA
        !

        call vocabu(NBOUN_INTE_1DIM,0_ip,0_ip)
        pari1 => kfl_fixbo_nsa
        call par_mygather()
        !
        ! BVNAT_NSA
        !
        call vocabu(NBOUN_REAL_2DIM,ndofn_nsa,0_ip)
        parr2 => bvnat_nsa 
        call par_mygather()
        !
        ! KFL_FIXRS_NSA
        !
        call vocabu(NBOPO_INTE_1DIM,0_ip,0_ip)
        pari1 => kfl_fixrs_nsa
        call par_mygather()
        !
        ! BVESS_NSA
        !
        call vocabu(NPOIN_REAL_2DIM,ndofn_nsa,0_ip)
        parr2 => bvess_nsa(1:ndofn_nsa,1:npoin,1)     
        call par_mygather()

        if(kfl_conbc_nsa==0) then
           call vocabu(NPOIN_REAL_2DIM,ndofn_nsa,0_ip)
           parr2 => bvess_nsa(1:ndofn_nsa,1:npoin,2)
           call par_mygather()          
           !
           ! KFL_FUNNO_NSA
           !
           call vocabu(NPOIN_INTE_1DIM,0_ip,0_ip)
           pari1 => kfl_funno_nsa        
           call par_mygather()
        end if

        !
        ! METEO reference arrays and initialization arrays:
        !
        
        ! q_nsa
        call vocabu(NPOIN_REAL_2DIM,ndofn_nsa,0_ip)
        parr2 => q_nsa(1:ndime+6,1:npoin)     
        call par_mygather()

        ! qref_nsa
        call vocabu(NPOIN_REAL_2DIM,ndofn_nsa,0_ip)
        parr2 => q_nsa(1:ndime+6,1:npoin)     
        call par_mygather()


        !
        ! REKEE_NSA
        !
        if (kfl_inkee_nsa(1) > 0 .or. kfl_inkee_nsa(2) > 0 .or. kfl_inkee_nsa(3) > 0 &
             .or. kfl_inkee_nsa(4) > 0 .or. kfl_inkee_nsa(5) > 0) then
           !        if(kfl_inkee_nsa==1) then

           if(kfl_benme_nsa >= 200) then
              call vocabu(NPOIN_REAL_2DIM,nkeep_nsa,0_ip)
              parr2 => rekee_nsa(1:nkeep_nsa+3,1:npoin)     
              call par_mygather()
           else
              call vocabu(NPOIN_REAL_2DIM,nkeep_nsa,0_ip)
              parr2 => rekee_nsa(1:nkeep_nsa,1:npoin)     
              call par_mygather()
           end if

        end if

        !
        ! Sponge (Mariano's obsolete version)
        !
        if(kfl_spong_nsa==1 .or. kfl_sponge_nsa > 0) then  
           call vocabu(NPOIN_REAL_2DIM,ndofn_nsa+1,0_ip)                                           
           parr2 => bspon_nsa(1:ndofn_nsa+1,1:npoin)                     
           call par_mygather()                                                                                          
        end if

        !
        ! Deallocate memory of Master
        !
        if(kfl_paral==0) call nsa_membcs(three)


     end if

  end select

  npari=0
  nparr=0
  nparc=0


end subroutine nsa_sendat
