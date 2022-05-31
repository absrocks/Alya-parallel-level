subroutine tur_memall()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_memall
  ! NAME 
  !    tur_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    turbulence equations      
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_turbul
  use mod_memory
  use def_kermod, only : kfl_adj_prob,kfl_ndvars_opt
  use mod_ADR, only : FULL_OSS
  use mod_ADR, only : A_OSS  
  use mod_ADR, only : AR_OSS 
  use mod_ADR, only : BUBBLE
  use mod_ADR, only : ADR_initialize_type
  use mod_ADR, only : ADR_check_and_compute_data
  use mod_ADR, only : ADR_allocate_projections_bubble_sgs 
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  implicit none
  integer(ip) :: iturb,ielem,pelty,pnode

  if( INOTMASTER ) then
     !
     ! Turbulence unknowns
     !
     call memory_alloca(mem_modul(1:2,modul),'UNTUR','tur_memall',untur,nturb_tur,npoin,ncomp_tur)
     call memory_alloca(mem_modul(1:2,modul),'TURMU','tur_memall',turmu,npoin)
     !
     ! DUNKN_TUR: Delta unknown for Aitken relaxation strategy
     !
     if( kfl_relax_tur == 2 ) then
        call memory_alloca(mem_modul(1:2,modul),'DUNKN_TUR','tur_memall',dunkn_tur,nturb_tur*npoin)        
     end if
     !
     ! UNOLD_TUR: Old solution
     !
     if(    output_postprocess_check_variable_postprocess(15_ip) .or.&
          & output_postprocess_check_variable_postprocess(16_ip) .or.&
          & output_postprocess_check_variable_postprocess(17_ip) ) then
        call memory_alloca(mem_modul(1:2,modul),'UNOLD_TUR','tur_memall',unold_tur,nturb_tur+1_ip,npoin)
     end if
     !
     ! UNPRO_TUR: Residual projections for OSS methods
     !
     if( kfl_ortho_tur >= 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'UNPRO_TUR','tur_memall',unpro_tur,nturb_tur,npoin)
        if( kfl_ortho_tur == 2 ) then  ! split oss
           call memory_alloca(mem_modul(1:2,modul),'UNPRR_TUR','tur_memall',unprr_tur,nturb_tur,npoin)
        end if
     end if
     !
     ! UNPGR_TUR: Residual projections for Shock capturing methos
     !
     if (kfl_shock_tur/=0)  call memory_alloca(mem_modul(1:2,modul),'UNPGR_TUR','tur_memall',unpgr_tur,nturb_tur,ndime, npoin)

     !
     ! DETUR_TUR, VITUR_TUR: projected variable density and viscosity
     !
     if( kfl_colev_tur >= 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'DETUR_TUR','tur_memall',detur_tur,npoin)
        call memory_alloca(mem_modul(1:2,modul),'VITUR_TUR','tur_memall',vitur_tur,npoin)
     end if
     !
     ! Solver memory
     !
     solve_sol => solve(1:1)
     call soldef(4_ip)
     solve_sol => solve(3:3)
     call soldef(4_ip)
     if(kfl_algor_tur==1) then
        solve_sol => solve(2:nturb_tur)
        call soldef(4_ip)
     end if
     !
     ! only for DDES, postprocessing value
     !
     if( kfl_ddesm_tur >= 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'FDDES_TUR','tur_memall',fddes_tur,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GDDES_TUR','tur_memall',gddes_tur,npoin)
     end if
     !
     ! only for SST, postprocessing values
     !
     if( TUR_SST_K_OMEGA) then
        call memory_alloca(mem_modul(1:2,modul),'SSTF1_TUR','tur_memall',sstf1_tur,npoin)
        call memory_alloca(mem_modul(1:2,modul),'SSTF2_TUR','tur_memall',sstf2_tur,npoin)
     end if
     if( TUR_SST_K_OMEGA .and. kfl_sasim_tur == 1) then
        call memory_alloca(mem_modul(1:2,modul),'SASSO_TUR','tur_memall',sasso_tur,npoin)
     end if
     if (TUR_FAMILY_K_EPS) &
          call memory_alloca(mem_modul(1:2,modul),'TUR_MAX_MIXLEN','tur_memall',tur_max_mixlen,npoin)
     !
     ! PRODU_TUR: production term
     !
     if( kfl_produ_tur == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'PRODU_TUR','tur_memall',produ_tur,npoin)
     end if
     !
     ! averaged values
     !
     if( output_postprocess_check_variable_postprocess(42_ip) ) then
        call memory_alloca(mem_modul(1:2,modul),'OLDED_TUR','tur_memall',olded_tur,npoin)
        call memory_alloca(mem_modul(1:2,modul),'AVTVI_TUR','tur_memall',avtvi_tur,npoin)
        olded_tur = 0.0_rp
        avtvi_tur = 0.0_rp
     end if
     if( output_postprocess_check_variable_postprocess(40_ip) ) then
        call memory_alloca(mem_modul(1:2,modul),'AVKEY_TUR','tur_memall',avkey_tur,npoin)
        avkey_tur = 0.0_rp
     end if
     if( output_postprocess_check_variable_postprocess(41_ip) ) then
        call memory_alloca(mem_modul(1:2,modul),'AVOME_TUR','tur_memall',avome_tur,npoin)
        avome_tur = 0.0_rp
     end if
     ! provisional variable (MATIAS)
     call memory_alloca(mem_modul(1:2,modul),'TURVI_TUR','tur_memall',turvi_tur,2_ip, mgaus, nelem)
     turvi_tur = 0.0_rp
     !
     ! adjoint materials
     !
     if(kfl_adj_prob == 1 ) then     
        call memory_alloca(mem_modul(1:2,modul),'UNTUR_FORW','tur_memall',untur_forw,nturb_tur,npoin,ncomp_tur)
        call memory_alloca(mem_modul(1:2,modul),'RESDIFF_TUR','tur_memall',resdiff_tur,kfl_ndvars_opt,npoin)

        if (nturb_tur == 1) then
           call memory_alloca(mem_modul(1:2,modul),'Rhsadjtur_tur','tur_memall',Rhsadjtur_tur,nelem)
           do ielem=1,nelem
	      pelty = ltype(ielem)
	      pnode = nnode(pelty)
	      call memory_alloca(mem_modul(1:2,modul),'Rhsadjtur_tur','tur_memall',Rhsadjtur_tur(ielem)%a,nturb_tur,pnode)
           end do
        endif

        call memory_alloca(mem_modul(1:2,modul),'RhsadjNas_tur','tur_memall',RhsadjNas_tur,nelem)

        do ielem=1,nelem
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           call memory_alloca(mem_modul(1:2,modul),'RhsadjNas_tur','tur_memall',RhsadjNas_tur(ielem)%a,ndime,pnode)
        end do

     endif !kfl_adj_prob

  else
     !
     ! Master allocate minimum memory
     !
     call memory_alloca(mem_modul(1:2,modul),'UNTUR','tur_memall',untur,4_ip,1_ip,3_ip)
     call memory_alloca(mem_modul(1:2,modul),'TURMU','tur_memall',turmu,1_ip)
     if( kfl_relax_tur == 2 ) then
        call memory_alloca(mem_modul(1:2,modul),'DUNKN_TUR','tur_memall',dunkn_tur,1_ip)      
     end if
     if (TUR_FAMILY_K_EPS) &
          call memory_alloca(mem_modul(1:2,modul),'TUR_MAX_MIXLEN','tur_memall',tur_max_mixlen,1_ip)
  end if
  !
  ! ADR type
  ! 
  do iturb = 1,nturb_tur
     call ADR_initialize_type(ADR_tur(iturb))
     ADR_tur(iturb) % kfl_time_integration   =  kfl_timei_tur
     ADR_tur(iturb) % kfl_time_step_strategy =  kfl_timco
     ADR_tur(iturb) % kfl_stabilization      =  kfl_ortho_tur
     ADR_tur(iturb) % kfl_shock              =  kfl_shock_tur
     ADR_tur(iturb) % kfl_time_lumped        =  0
     ADR_tur(iturb) % kfl_tau_strategy       =  kfl_taust_tur
     ADR_tur(iturb) % kfl_laplacian          =  0 
     ADR_tur(iturb) % kfl_nonlinear_sgs      =  kfl_sgsno_tur 
     ADR_tur(iturb) % kfl_time_sgs           =  kfl_sgsti_tur
     ADR_tur(iturb) % kfl_time_bubble        =  kfl_tibub_tur
     ADR_tur(iturb) % kfl_time_scheme        =  kfl_tisch_tur
     ADR_tur(iturb) % kfl_time_order         =  kfl_sgsac_tur
     ADR_tur(iturb) % kfl_manufactured       =  kfl_exacs_tur
     ADR_tur(iturb) % kfl_length             =  kfl_ellen_tur
     ADR_tur(iturb) % number_euler_steps     =  neule_tur

     ADR_tur(iturb) % lun_output4            =  int(momod(modul) % lun_outpu,4)
     ADR_tur(iturb) % bemol                  =  bemol_tur
     ADR_tur(iturb) % tau_parameters(1:3)    =  staco_tur(1:3)
     ADR_tur(iturb) % shock                  =  shock_tur  
     call ADR_check_and_compute_data(ADR_tur(iturb))
     call ADR_allocate_projections_bubble_sgs(ADR_tur(iturb))
  end do

end subroutine tur_memall
 
