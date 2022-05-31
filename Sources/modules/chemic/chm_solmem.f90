subroutine chm_solmem()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_memall
  ! NAME 
  !    chm_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master 
  use def_domain
  use def_solver
  use def_chemic
  use mod_memory
  use mod_ADR,             only : FULL_OSS
  use mod_ADR,             only : A_OSS  
  use mod_ADR,             only : AR_OSS 
  use mod_ADR,             only : BUBBLE
  use mod_ADR,             only : ADR_initialize_type
  use mod_ADR,             only : ADR_check_and_compute_data
  use mod_ADR,             only : ADR_allocate_projections_bubble_sgs
  use mod_chm_rk_explicit, only : chm_rk_explicit_memory
  use mod_chm_sectional_soot_model, only : chm_sectional_soot_memory_ssm 
  use mod_chm_sectional_soot_model, only : chm_initialization_ssm 

  implicit none
  integer(ip) :: ielem,pelty,pgaus, ipoin,ireac

  if( INOTEMPTY ) then

     call memory_alloca(mem_modul(1:2,modul),'SPHEC',    'chm_solmem',sphec,npoin,6_ip,2_ip)
     call memory_alloca(mem_modul(1:2,modul),'AVDEN_CHM','chm_solmem',avden_chm,npoin)
     call memory_alloca(mem_modul(1:2,modul),'AVMSK_CHM','chm_solmem',avmsk_chm,npoin)        ! Average mass source from spray
     call memory_alloca(mem_modul(1:2,modul),'MASSK',    'chm_solmem',massk,npoin,nclas_chm)

     !
     ! Kernel shared variables
     !
     call memory_alloca(mem_modul(1:2,modul),'WMEAN','chm_solmem',wmean,npoin,ncomp_chm)

     !
     ! Projection of dt/rho 
     !
     call memory_alloca(mem_modul(1:2,modul),'DT_RHO_CHM','chm_solmem',dt_rho_chm,npoin)
     if (kfl_spray_chm /= 0_ip ) call memory_alloca(mem_modul(1:2,modul),'DT_CHM','chm_solmem',dt_chm,npoin)

     call memory_alloca(mem_modul(1:2,modul),'CHEMICAL_HEAT','chm_solmem',chemical_heat,nelem)
     do ielem = 1,nelem
        pelty = ltype(ielem)
        pgaus = ngaus(pelty)
        call memory_alloca(mem_modul(1:2,modul),'CHEMICAL_HEAT % A','chm_solmem',chemical_heat(ielem)%a,pgaus,1_ip,1_ip)
     end do

     call memory_alloca(mem_modul(1:2,modul),'RADIATIVE_HEAT','chm_solmem',radiative_heat,nelem)
     do ielem = 1,nelem 
        pelty = ltype(ielem)  
        pgaus = ngaus(pelty)
        call memory_alloca(mem_modul(1:2,modul),'RADIATIVE_HEAT % A','chm_solmem',radiative_heat(ielem)%a,pgaus,1_ip,1_ip)
     end do

     !
     ! Allocation transport properties
     !  
     if (kfl_lookg_chm > 0) then 

        if (kfl_model_chm /= 4) then

           call memory_alloca(mem_modul(1:2,modul),'CONDU_GP','chm_solmem',condu_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'CONDU_GP','chm_solmem',condu_gp(ielem)%a,pgaus,1_ip,1_ip)
              condu_gp(ielem)%a=0.0_rp
           end do

           call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP','chm_solmem',sphec_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP','chm_solmem',sphec_gp(ielem)%a,pgaus,1_ip,1_ip)
              sphec_gp(ielem)%a=1.0_rp
           end do

           call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp(ielem)%a,pgaus,1_ip,1_ip)
              visco_gp(ielem)%a=0.0_rp
           end do

           call memory_alloca(mem_modul(1:2,modul),'WMEAN_GP','chm_solmem',wmean_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'WMEAN_GP','chm_solmem',wmean_gp(ielem)%a,pgaus,4_ip,1_ip)
              wmean_gp(ielem)%a=0.0_rp
           end do

           call memory_alloca(mem_modul(1:2,modul),'TEMPE_GP','chm_solmem',tempe_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'TEMPE_GP','chm_solmem',tempe_gp(ielem)%a,pgaus,1_ip,1_ip)
              tempe_gp(ielem)%a=200.0_rp
           end do

           call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_HT','chm_solmem',sphec_gp_ht,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_HT','chm_solmem',sphec_gp_ht(ielem)%a,pgaus,6_ip,1_ip)
              sphec_gp_ht(ielem)%a=1.0_rp
           end do

           call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_LT','chm_solmem',sphec_gp_lt,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_LT','chm_solmem',sphec_gp_lt(ielem)%a,pgaus,6_ip,1_ip)
              sphec_gp_lt(ielem)%a=1.0_rp
           end do

           call memory_alloca(mem_modul(1:2,modul),'MASS_GP','chm_solmem',mass_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'MASS_GP','chm_solmem',mass_gp(ielem)%a,pgaus,nclas_chm,1_ip)
              mass_gp(ielem)%a=0.0_rp
           end do

           call memory_alloca(mem_modul(1:2,modul),'RRT_GP','chm_solmem',rrt_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'RRT_GP','chm_solmem',rrt_gp(ielem)%a,pgaus,nclas_chm,1_ip)
              rrt_gp(ielem)%a=0.0_rp
           end do

        end if

        if ( kfl_model_chm == 3 ) then
           call memory_alloca(mem_modul(1:2,modul),'HK_GP','chm_solmem',hk_gp,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'HK_GP','chm_solmem',hk_gp(ielem)%a,nclas_chm,pgaus,1_ip)
              hk_gp(ielem)%a=0.0_rp
           end do
        end if

     else
        call memory_alloca(mem_modul(1:2,modul),'VISCK','chm_solmem',visck,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'CONDK','chm_solmem',condk,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'SPHEK','chm_solmem',sphek,npoin,nclas_chm)
     endif

     !
     ! Finite rate chemistry model variables
     !
     if ( kfl_model_chm == 3 ) then
        call memory_alloca(mem_modul(1:2,modul),'ENTHALPY_TRANSPORT','chm_solmem',enthalpy_transport,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'ENTHALPY_TRANSPORT','chm_solmem',enthalpy_transport(ielem)%a,pgaus,ndime,1_ip)
           enthalpy_transport(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'DIV_ENTHALPY_TRANSPORT','chm_solmem',div_enthalpy_transport,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'DIV_ENTHALPY_TRANSPORT','chm_solmem',div_enthalpy_transport(ielem)%a,pgaus,1_ip,1_ip)
           div_enthalpy_transport(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'dummy_enthalpy','chm_solmem',dummy_enthalpy,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'dummy_enthalpy','chm_solmem',dummy_enthalpy(ielem)%a,pgaus,1_ip,1_ip)
           dummy_enthalpy(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'GRAD_YK',                 'chm_solmem',grad_Yk,                 nspec_chm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_YK_HK',              'chm_solmem',grad_Yk_hk,              nspec_chm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_T',                  'chm_solmem',grad_T,                  ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AUX_NODES',               'chm_solmem',aux_nodes,               npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_H',                  'chm_solmem',grad_H,                  ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_K_SINGLE',           'chm_solmem',grad_k_single,           ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'ENTHALPY_TRANSPORT_NODES','chm_solmem',enthalpy_transport_nodes,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'Y_K_N',                   'chm_solmem',Y_k_n,                   npoin,nspec_chm,ncomp_chm)
        call memory_alloca(mem_modul(1:2,modul),'REACT_IND',               'chm_solmem',React_ind,               npoin,nreac_chm)
        call memory_alloca(mem_modul(1:2,modul),'CORR_CHM',                'chm_solmem',Corr_chm,                npoin)
        call memory_alloca(mem_modul(1:2,modul),'SRC_CHM',                 'chm_solmem',src_chm,                 npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'ELEM_C',                  'chm_solmem',elem_c,                  npoin)
        call memory_alloca(mem_modul(1:2,modul),'ELEM_H',                  'chm_solmem',elem_h,                  npoin)
        call memory_alloca(mem_modul(1:2,modul),'ELEM_N',                  'chm_solmem',elem_n,                  npoin)
        call memory_alloca(mem_modul(1:2,modul),'ELEM_O',                  'chm_solmem',elem_o,                  npoin)
        call memory_alloca(mem_modul(1:2,modul),'MIXFR_CHM',               'chm_solmem',mixfr_chm,               npoin)
        call memory_alloca(mem_modul(1:2,modul),'PROG_VAR_CHM',            'chm_solmem',prog_var_chm,            npoin)
        call memory_alloca(mem_modul(1:2,modul),'SUM_REAC_CHM',            'chm_solmem',sum_reac_chm,            npoin)
        call memory_alloca(mem_modul(1:2,modul),'HRR_CHM',                 'chm_solmem',hrr_chm,                 npoin)
        call memory_alloca(mem_modul(1:2,modul),'HRR_AVG_CHM',             'chm_solmem',hrr_avg_chm,             npoin)
        call memory_alloca(mem_modul(1:2,modul),'ENTHA_CHM',               'chm_solmem',entha_chm,               nclas_chm,npoin)

        do ipoin = 1,npoin
           do ireac = 1,nreac_chm
              React_ind(ipoin,ireac)       = 1_ip
           enddo
           Corr_chm(ipoin)        = 1_ip
        enddo


     elseif( kfl_model_chm == 2 .or. kfl_model_chm == 1) then
        !
        ! Species fields for radiation for flamelet model
        !
        call memory_alloca(mem_modul(1:2,modul),'RSPEC_CHM'   ,'chm_solmem',rspec_chm,2_ip,npoin)   
        call memory_alloca(mem_modul(1:2,modul),'ZGRADMAX_CHM','chm_solmem',zgradmax_chm,npoin)     
        call memory_alloca(mem_modul(1:2,modul),'PHI_CHM'     ,'chm_solmem',phi_chm,npoin)          
        call memory_alloca(mem_modul(1:2,modul),'AVCHM_CHM'   ,'chm_solmem',avchm_chm,npoin)        

        call memory_alloca(mem_modul(1:2,modul),'LESCL','chm_solmem',lescl,npoin)

        call memory_alloca(mem_modul(1:2,modul),'AVY_CHM'  ,'chm_solmem',avY_chm, npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVZ_CHM'  ,'chm_solmem',avZ_chm, npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVYV_CHM' ,'chm_solmem',avYV_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVZV_CHM' ,'chm_solmem',avZV_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVZ2_CHM' ,'chm_solmem',avZ2_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVY2_CHM' ,'chm_solmem',avY2_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVL_CHM'  ,'chm_solmem',avL_chm, npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVL2_CHM' ,'chm_solmem',avL2_chm,npoin)

        !
        ! Allocation scalar dissipation rates
        !  
        call memory_alloca(mem_modul(1:2,modul),'XYR_CHM','chm_solmem',xYr_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'XZR_CHM','chm_solmem',xZr_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'XYS_CHM','chm_solmem',xYs_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'XZS_CHM','chm_solmem',xZs_chm,npoin)

        call memory_alloca(mem_modul(1:2,modul),'AVXYR_CHM',    'chm_solmem',avxYr_chm,     npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVXZR_CHM',    'chm_solmem',avxZr_chm,     npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVXYS_CHM',    'chm_solmem',avxYs_chm,     npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVXZS_CHM',    'chm_solmem',avxZs_chm,     npoin)

        if (associated(posttable_fw)) then
           call memory_alloca(mem_modul(1:2,modul),'AVPOSTTAB_CHM','chm_solmem',avposttab_chm, npoin, posttable_fw % main_table % nvar)
        endif

     endif


     ! Allocate variables for CMC model
     if ( kfl_model_chm == 4 ) then
        print*, '--| ALYA     CHEMIC: ALLOCATING MEMORY FOR CMC MODEL...'
        nullify(enthalp_CMC_chm)
        nullify(temp_CMC_chm)
        nullify(Yk_CMC_chm)
        nullify(src_Yk_CMC_chm)
        nullify(Yk_int_CMC_chm)
        nullify(densi_int_CMC_chm)
        nullify(visco_lam_int_CMC_chm)
        nullify(enthalp_int_CMC_chm)
        nullify(temp_int_CMC_chm)
        nullify(src_Yk_int_CMC_chm)
        nullify(hrr_CMC_chm)
        nullify(deriv2_Yk_CMC_chm)
        nullify(veloc_CFD_chm)
        nullify(Zavg_CFD_chm)
        nullify(Zvar_CFD_chm)
        nullify(Xtot_CFD_chm)
        nullify(grad_Zavg_CFD_chm)
        nullify(visco_turb_CFD_chm)
        nullify(kfl_bc_type_spec_CMC_chm)
        nullify(W_k)
        nullify(coeff_cp_k)
        nullify(hrr_chm)
        nullify(hrr_avg_chm)
        nullify(condk)
        nullify(sphek)

        if (kfl_solve_enth_CMC_chm == 0) then
           call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_ENTHALPY_CMC_CHM','chm_solmem',enthalp_CMC_chm,nZ_CMC_chm,1_ip)
        else
           call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_ENTHALPY_CMC_CHM','chm_solmem',enthalp_CMC_chm,nZ_CMC_chm,npoin)
        end if
        call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_TEMPERATURE_CMC_CHM','chm_solmem',temp_CMC_chm,nZ_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_MASS_FRACTIONS_CMC_CHM','chm_solmem',Yk_CMC_chm,nZ_CMC_chm,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'CONDITIONAL_OMEGA_MASS_FRACT_CMC_CHM','chm_solmem',src_Yk_CMC_chm,nZ_CMC_chm,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_MASS_FRACTIONS_CMC_CHM','chm_solmem',Yk_int_CMC_chm,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_DENSITY_CMC_CHM','chm_solmem',densi_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_LAMINAR_VISC_CMC_CHM','chm_solmem',visco_lam_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_ENTHALPY_CMC_CHM','chm_solmem',enthalp_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_TEMPERATURE_CMC_CHM','chm_solmem',temp_int_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_OMEGA_MASS_FRACT_CMC_CHM','chm_solmem',src_Yk_int_CMC_chm,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'HEAT_RELEASE_CMC_CHM','chm_solmem',hrr_CMC_chm,nZ_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'DERIV2_Yk_CMC_CHM','chm_solmem',deriv2_Yk_CMC_chm,nZ_CMC_chm,npoin,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'VELOCITY_CFD_CHM','chm_solmem',veloc_CFD_chm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'MIXT_FRAC_CFD_CHM','chm_solmem',Zavg_CFD_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'MIXT_FRAC_VAR_CFD_CHM','chm_solmem',Zvar_CFD_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'TOTAL_SCAL_DISSIP_RATE_CFD_CHM','chm_solmem',Xtot_CFD_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'MIXTURE_FRACTION_GRADIENT_CFD_CHM','chm_solmem',grad_Zavg_CFD_chm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'MASS_TURBULENT_DIFFUSION_COEFF_CFD_CHM','chm_solmem',visco_turb_CFD_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'TYPE_BC_SPECIES_CMC_CHM','chm_solmem',kfl_bc_type_spec_CMC_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'W_K','chm_solmem',W_k,nclas_chm)
        call memory_alloca(mem_modul(1:2,modul),'COEFF_CP_K','chm_solmem',coeff_cp_k,nclas_chm,15_ip)
        call memory_alloca(mem_modul(1:2,modul),'HRR_CHM','chm_solmem',hrr_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'HRR_AVG_CHM','chm_solmem',hrr_avg_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_CONDUCTIVITY_CHM','chm_solmem',condk,npoin,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'UNCONDITIONAL_SPECIFIC_HEAT_CHM','chm_solmem',sphek,npoin,1_ip)



        ! Initialize to 0 previous matrices
        if (kfl_solve_enth_CMC_chm == 0) then
           enthalp_CMC_chm(1:nZ_CMC_chm,1)                   = 0.0_rp
        else
           enthalp_CMC_chm(1:nZ_CMC_chm,1:npoin)             = 0.0_rp
        end if
        temp_CMC_chm(1:nZ_CMC_chm,1:npoin)                   = 300.0_rp
        Yk_CMC_chm(1:nZ_CMC_chm,1:npoin,1:nclas_chm)         = 0.0_rp
        src_Yk_CMC_chm(1:nZ_CMC_chm,1:npoin,1:nclas_chm)     = 0.0_rp
        Yk_int_CMC_chm(1:npoin,1:nclas_chm)                  = 0.0_rp
        densi_int_CMC_chm(1:npoin)                           = 0.0_rp
        visco_lam_int_CMC_chm(1:npoin)                       = 0.0_rp
        enthalp_int_CMC_chm(1:npoin)                         = 0.0_rp
        temp_int_CMC_chm(1:npoin)                            = 0.0_rp
        src_Yk_int_CMC_chm(1:npoin,1:nclas_chm)              = 0.0_rp
        hrr_CMC_chm(1:nZ_CMC_chm,1:npoin)                    = 0.0_rp
        deriv2_Yk_CMC_chm(1:nZ_CMC_chm,1:npoin,1:nclas_chm)  = 0.0_rp
        !!!!!!!!!!!!!!!! PROVISIONAL
        veloc_CFD_chm(1:ndime,1:npoin)                       = 0.0_rp
        Zavg_CFD_chm(1:npoin)                                = 0.5_rp
        Zvar_CFD_chm(1:npoin)                                = 0.05_rp
        !!!!Xtot_CFD_chm(1:npoin)                                = 20.0_rp
        Xtot_CFD_chm(1:npoin)                                = 0.0_rp
        grad_Zavg_CFD_chm(1:ndime,1:npoin)                   = 0.0_rp
        visco_turb_CFD_chm(1:npoin)                          = 0.0_rp
        !!!!!!!!!!!!!!! PROVIISONAL
        kfl_bc_type_spec_CMC_chm(1:npoin)                    = 1
        W_k(1:nclas_chm)                                     = 0.0_rp
        coeff_cp_k(1:nclas_chm,1:15)                         = 0.0_rp
        hrr_chm(1:npoin)                                     = 0.0_rp
        hrr_avg_chm(1:npoin)                                 = 0.0_rp
        condk(1:npoin,1)                                     = 0.0_rp
        sphek(1:npoin,1)                                     = 0.0_rp


        if (.not. associated(therm)) then
           nullify(therm)
           call memory_alloca(mem_modul(1:2,modul),'THERM_CMC_CHM','chm_solmem',therm,npoin,ncomp_chm)
           therm(1:npoin,1:ncomp_chm) = 0.0_rp
        end if

        if (.not. associated(tempe)) then
           nullify(tempe)
           call memory_alloca(mem_modul(1:2,modul),'THERM_CMC_CHM','chm_solmem',tempe,npoin,ncomp_chm)
           tempe(1:npoin,1:ncomp_chm) = 0.0_rp
        end if


        if (kfl_solve_enth_CMC_chm /= 0) then
           nullify(deriv2_enthalp_CMC_chm)
           call memory_alloca(mem_modul(1:2,modul),'DERIV2_ENTHALP_CMC_CHM','chm_solmem',deriv2_enthalp_CMC_chm,nZ_CMC_chm,npoin)
           deriv2_enthalp_CMC_chm(1:nZ_CMC_chm,1:npoin) = 0.0_rp
        end if


        call memory_alloca(mem_modul(1:2,modul),'CONDU_GP_CMC_CHM','chm_solmem',condu_gp_CMC_chm,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'CONDU_GP_CMC_CHM','chm_solmem',condu_gp_CMC_chm(ielem)%a,pgaus,1_ip,1_ip)
           condu_gp_CMC_chm(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_CMC_CHM','chm_solmem',sphec_gp_CMC_chm,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_CMC_CHM','chm_solmem',sphec_gp_CMC_chm(ielem)%a,pgaus,1_ip,1_ip)
           sphec_gp_CMC_chm(ielem)%a=1.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'VISCO_GP_CMC_CHM','chm_solmem',visco_gp_CMC_chm,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'VISCO_GP_CMC_CHM','chm_solmem',visco_gp_CMC_chm(ielem)%a,nZ_CMC_chm,pgaus,1_ip)
           visco_gp_CMC_chm(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp(ielem)%a,pgaus,1_ip,1_ip)
           visco_gp(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'SPEC_VOL_GP_CMC_CHM','chm_solmem',spvol_gp_CMC_chm,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'SPEC_VOL_GP_CMC_CHM','chm_solmem',spvol_gp_CMC_chm(ielem)%a,nZ_CMC_chm,pgaus,1_ip)
           spvol_gp_CMC_chm(ielem)%a=0.0_rp
        end do

        call memory_alloca(mem_modul(1:2,modul),'DENSI_GP','chm_solmem',densi_gp,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'DENSI_GP','chm_solmem',densi_gp(ielem)%a,pgaus,1_ip,1_ip)
           densi_gp(ielem)%a=0.0_rp
        end do

     end if


     !
     ! Allocate spray terms
     !
     if ( kfl_spray_chm /= 0 ) then

        call memory_alloca(mem_modul(1:2,modul),'AVS_CHM'  ,'chm_solmem',avS_chm  ,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVS0_CHM' ,'chm_solmem',avS0_chm ,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVD32_CHM','chm_solmem',avd32_chm,npoin)

        call memory_alloca(mem_modul(1:2,modul),'SIGMA_CHM','chm_solmem',Sigma_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'SIGM0_CHM','chm_solmem',Sigm0_chm,npoin)
        call memory_alloca(mem_modul(1:2,modul),'D32_CHM'  ,'chm_solmem',d32_chm,  npoin)

        call memory_alloca(mem_modul(1:2,modul),'SIGMA_GP_CHM' ,'chm_solmem',sigma_gp_chm ,nelem)
        call memory_alloca(mem_modul(1:2,modul),'SIGMA0_GP_CHM','chm_solmem',sigma0_gp_chm,nelem)
        call memory_alloca(mem_modul(1:2,modul),'D32_GP_CHM'   ,'chm_solmem',d32_gp_chm   ,nelem)

        do ielem = 1,nelem 
           pelty = ltype(ielem)  
           pgaus = ngaus(pelty)      
           call memory_alloca(mem_modul(1:2,modul),'SIGMA_GP_CHM(IELEM)' ,'chm_solmem',sigma_gp_chm(ielem)%a ,pgaus,1_ip,1_ip)
           call memory_alloca(mem_modul(1:2,modul),'SIGMA0_GP_CHM(IELEM)','chm_solmem',sigma0_gp_chm(ielem)%a,pgaus,1_ip,1_ip)
           call memory_alloca(mem_modul(1:2,modul),'D32_GP_CHM(IELEM)'   ,'chm_solmem',d32_gp_chm(ielem)%a   ,pgaus,1_ip,1_ip)
        end do

     end if
     !
     ! Allocate variables for level set
     !
     if ( kfl_spray_chm == 2 ) then
        call memory_alloca(mem_modul(1:2,modul),'LAP_PHI_LEVSET_CHM' ,'chm_solmem',lap_phi_levSet_chm ,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_PHI_LEVSET_CHM','chm_solmem',grad_phi_levSet_chm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_PHIC_LEVSET_CHM','chm_solmem',grad_phic_levSet_chm,ndime,npoin)
     end if

     !! DMM !
     !! DMM ! Allocate variables for soot model
     !! DMM !
     !! DMM if ( kfl_soot_chm /= 0 ) then 
     !! DMM    call chm_sectional_soot_memory_ssm()
     !! DMM    call chm_initialization_ssm()
     !! DMM endif

  else
     !
     ! Projection of dt/rho 
     !
     call memory_alloca(mem_modul(1:2,modul),'DT_RHO_CHM','chm_solmem',dt_rho_chm,1_ip)
     if (kfl_spray_chm /= 0_ip ) call memory_alloca(mem_modul(1:2,modul),'DT_CHM','chm_solmem',dt_chm,1_ip)

     !
     ! Allocation scalar dissipation rates
     !  
     call memory_alloca(mem_modul(1:2,modul),'XYR_CHM','chm_solmem',xYr_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'XZR_CHM','chm_solmem',xZr_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'XYS_CHM','chm_solmem',xYs_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'XZS_CHM','chm_solmem',xZs_chm,1_ip)

     call memory_alloca(mem_modul(1:2,modul),'AVXYR_CHM','chm_solmem',avxYr_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'AVXZR_CHM','chm_solmem',avxZr_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'AVXYS_CHM','chm_solmem',avxYs_chm,1_ip)
     call memory_alloca(mem_modul(1:2,modul),'AVXZS_CHM','chm_solmem',avxZs_chm,1_ip)


     if (kfl_lookg_chm > 0) then 
        call memory_alloca(mem_modul(1:2,modul),'CONDU_GP','chm_solmem',condu_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP','chm_solmem',sphec_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'WMEAN_GP','chm_solmem',wmean_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'TEMPE_GP','chm_solmem',tempe_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_HT','chm_solmem',sphec_gp_ht,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_LT','chm_solmem',sphec_gp_lt,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'MASS_GP','chm_solmem',mass_gp,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'RRT_GP','chm_solmem',rrt_gp,1_ip)

        call memory_alloca(mem_modul(1:2,modul),'CONDU_GP','chm_solmem',condu_gp(1)%a,1_ip,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP','chm_solmem',sphec_gp(1)%a,1_ip,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'VISCO_GP','chm_solmem',visco_gp(1)%a,1_ip,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'WMEAN_GP','chm_solmem',wmean_gp(1)%a,1_ip,4_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'TEMPE_GP','chm_solmem',tempe_gp(1)%a,1_ip,4_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_HT','chm_solmem',sphec_gp_ht(1)%a,1_ip,6_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SPHEC_GP_LT','chm_solmem',sphec_gp_lt(1)%a,1_ip,6_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'MASS_GP','chm_solmem',mass_gp(1)%a,1_ip,nclas_chm,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'RRT_GP','chm_solmem',rrt_gp(1)%a,1_ip,nclas_chm,1_ip)
     endif

     !
     ! Allocate spray terms
     !
     if ( kfl_spray_chm /= 0 ) then

        call memory_alloca(mem_modul(1:2,modul),'AVS_CHM'  ,'chm_solmem',avS_chm  ,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'AVS0_CHM' ,'chm_solmem',avS0_chm ,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'AVD32_CHM','chm_solmem',avd32_chm,1_ip)

        call memory_alloca(mem_modul(1:2,modul),'SIGMA_GP_CHM' ,'chm_solmem',sigma_gp_chm ,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SIGMA0_GP_CHM','chm_solmem',sigma0_gp_chm,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'D32_GP_CHM'   ,'chm_solmem',d32_gp_chm   ,1_ip)

     end if

     !
     ! Allocate variables for level set
     !
     if ( kfl_spray_chm == 2 ) then
        call memory_alloca(mem_modul(1:2,modul),'LAP_PHI_LEVSET_CHM' ,'chm_solmem',lap_phi_levSet_chm   ,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_PHI_LEVSET_CHM','chm_solmem',grad_phi_levSet_chm  ,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'GRAD_PHIC_LEVSET_CHM','chm_solmem',grad_phic_levSet_chm,1_ip,1_ip)
     end if

  end if

  if ( kfl_model_chm == 3 ) then
     call memory_alloca(mem_modul(1:2,modul),'COEFF_CP_K','chm_solmem',coeff_cp_k,nspec_chm,15_ip)
     call memory_alloca(mem_modul(1:2,modul),'W_K','chm_solmem',W_k,nspec_chm)
  endif

  !
  ! Allocate variables for soot model
  !
  if ( kfl_soot_chm /= 0 ) then 
     call chm_sectional_soot_memory_ssm()
     !!call chm_initialization_ssm()
  endif

  !
  ! Memory (Solver)
  !
  if ( kfl_model_chm == 4 ) then
     solve(1) % ndofn = nvar_CMC_chm
     solve(2) % ndofn = nvar_CMC_chm
  else
     solve(1) % ndofn = nclas_chm
     solve(2) % ndofn = nclas_chm
  end if

  solve_sol => solve(1:)
  call soldef(4_ip)
  !
  ! Residual
  !
  if (kfl_model_chm == 4) then
     call memory_alloca(mem_modul(1:2,modul),'RIPTS_CHM','chm_solmem',ripts_chm,nvar_CMC_chm*nclas_chm)
  else
     call memory_alloca(mem_modul(1:2,modul),'RIPTS_CHM','chm_solmem',ripts_chm,nclas_chm*nclas_chm)
  end if

  !
  ! Boundary conditions
  !
  solve(1) % bvess     => bvess_chm
  solve(1) % kfl_fixno => kfl_fixno_chm 
  solve(2) % bvess     => bvess_chm
  solve(2) % kfl_fixno => kfl_fixno_chm
  !  
  ! Memory
  !
  call chm_rk_explicit_memory()     

end subroutine chm_solmem

