subroutine tem_memall()
  !-----------------------------------------------------------------------
  !****f* temper/tem_memall
  ! NAME 
  !    tem_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    temperature equation
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_kermod, only : kfl_adj_prob, kfl_ndvars_opt
  use def_domain
  use def_solver
  use def_temper
  use mod_memchk
  use mod_memory
  use mod_ADR, only : FULL_OSS
  use mod_ADR, only : A_OSS  
  use mod_ADR, only : AR_OSS 
  use mod_ADR, only : BUBBLE
  use mod_ADR, only : ADR_initialize_type
  use mod_ADR, only : ADR_check_and_compute_data
  use mod_ADR, only : ADR_allocate_projections_bubble_sgs
  implicit none
  integer(ip) :: ielem,pelty,pnode
  integer(4)  :: istat, ncomp
  !
  ! Problem unknowns TEMPER, TESGS and solver initialization
  !
  if( INOTMASTER ) then
     !
     ! TEMPER: Temperature unknown 
     !
     call memory_alloca(mem_modul(1:2,modul),'TEMPE','tem_memall',tempe,max(1_ip,nunkn_tem),ncomp_tem)
     call memory_alloca(mem_modul(1:2,modul),'THERM','tem_memall',therm,max(1_ip,nunkn_tem),ncomp_tem)

     if (kfl_regim_tem /= 4) tempe => therm 
        !
        ! GRTEM: Temperature gradients
        !
     if(kfl_ellen_tem==-1) then
        call memory_alloca(mem_modul(1:2,modul),'GRTEM_TEM','tem_memall',grtem_tem,ndime,nunkn_tem)
        end if
        !
        ! TEOLD_TEM: Old temperaure
        !
        if(postp(1) % npp_stepi(8)>0) then
           call memory_alloca(mem_modul(1:2,modul),'TEOLD_TEM','tem_memall',teold_tem,nunkn_tem) 
        end if
        !
        ! AVTEM_TEM: Average temperature
        !
        if(postp(1) % npp_stepi(5)>0) then
           call memory_alloca(mem_modul(1:2,modul),'AVTEM_TEM','tem_memall',avtem_tem,nunkn_tem)
        end if
        !
        ! AVTE2_TEM: Average tempe**2
        !
        if(postp(1) % npp_stepi(18)>0) then
           call memory_alloca(mem_modul(1:2,modul),'AVTE2_TEM','tem_memall',avte2_tem,nunkn_tem)
        end if
        !
        ! AVTEV_TEM: Average tempe*veloc
        !
        if(postp(1) % npp_stepi(19)>0) then
           call memory_alloca(mem_modul(1:2,modul),'AVTEV_TEM','tem_memall',avtev_tem,ndime, nunkn_tem)
        end if
        !
        ! AVDEN_TEM: Average density
        !
        if(postp(1) % npp_stepi(20)>0) then
           call memory_alloca(mem_modul(1:2,modul),'AVDEN_TEM','tem_memall',avden_tem,nunkn_tem)
        end if
        !
        ! FVVEL_TEM: Average rho*veloc
        !
        if(postp(1) % npp_stepi(21)>0) then
           call memory_alloca(mem_modul(1:2,modul),'FVVEL_TEM','tem_memall',fvvel_tem,ndime, nunkn_tem)
        end if
        !
        ! AVRES_TEM: Average residual heat flux
        !
        if(postp(1) % npp_stepi(32)>0) then
           call memory_alloca(mem_modul(1:2,modul),'AVRES_TEM','tem_memall',avres_tem,nunkn_tem)
        end if
        !
        ! VELOC, TURMU: dynamic coupling
        !
        if(kfl_inter_tem==1) then
           if(kfl_advec_tem==1) then
              allocate(veloc(ndime,npoin,1),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'VELOC','tem_memall',veloc)
           end if
           if(kfl_cotur_tem==1) then
              allocate(turmu(npoin),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'TURMU','tem_memall',turmu)
           end if
        end if
        !
        ! Water vapor concentration gradients
        !
        call memory_alloca(mem_modul(1:2,modul),'GRADC_TEM','tem_memall',gradc_tem,ndime,nunkn_tem)
        !
        ! variables needed for adjoint
        !     
        if (kfl_adj_prob == 1) then
           !
           ! TEMPER_FORWARD: Temperature known
           !
           call memory_alloca(mem_modul(1:2,modul),'TEMPE_FORW','tem_memall',tempe_forw, nunkn_tem,ncomp_tem)
           !
           ! resdiff_tem: partial derivatives of R w.r.t design variables 
           !
           call memory_alloca(mem_modul(1:2,modul),'RESDIFF_TEM','tem_memall',resdiff_tem, kfl_ndvars_opt,solve_sol(1) % nzrhs)
           ! Coupling with nastin for coupled adjoint solution
           if(kfl_coupl(ID_TEMPER,ID_NASTIN) == 1 ) then
              call memory_alloca(mem_modul(1:2,modul),'RhsadjNas_tem','tem_memall',RhsadjNas_tem, nelem)
              do ielem=1,nelem
                pelty = ltype(ielem)
                pnode = nnode(pelty)
                call memory_alloca(mem_modul(1:2,modul),'RhsadjNas_tem','tem_memall',RhsadjNas_tem(ielem)%a, ndime,pnode)
              end do
           end if
           !
           ! Coupling with chemic for coupled adjoint solution
           !
           if(kfl_coupl(ID_TEMPER,ID_CHEMIC) == 1 ) then
             call memory_alloca(mem_modul(1:2,modul),'RhsadjChe_tem','tem_memall', RhsadjChe_tem, nelem)
             do ielem=1,nelem
                pelty = ltype(ielem)
                pnode = nnode(pelty)
                call memory_alloca(mem_modul(1:2,modul),'RhsadjChe_tem','tem_memall',RhsadjChe_tem(ielem)%a, pnode,nspec)
             end do
           endif
  
        endif
  else

     allocate(tempe(1,ncomp_tem),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'TEMPE','tem_memall',tempe)
     allocate(therm(1,ncomp_tem),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'THERM','tem_memall',therm)

     if (kfl_regim_tem /= 4) tempe => therm

     if(postp(1) % npp_stepi(5)>0) then
        allocate(avtem_tem(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'AVTEM_TEM','tem_memall',avtem_tem)        
     end if
     if(postp(1) % npp_stepi(32)>0) then
        allocate(avres_tem(1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'AVRES_TEM','tem_memall',avres_tem)
     end if
     if(kfl_ellen_tem==-1) then
        allocate(grtem_tem(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GRTEM_TEM','tem_memall',grtem_tem)        
     end if

  end if
  !
  ! Solver memory
  !
  solve_sol => solve(1:)
  call soldef(4_ip)
  !
  ! ADR type
  !
  call ADR_initialize_type(ADR_tem)
  ADR_tem % kfl_time_integration   =  kfl_timei_tem
  ADR_tem % kfl_time_step_strategy =  kfl_timco
  ADR_tem % kfl_stabilization      =  kfl_ortho_tem
  ADR_tem % kfl_shock              =  kfl_shock_tem
  ADR_tem % kfl_time_lumped        =  0
  ADR_tem % kfl_tau_strategy       =  kfl_taust_tem
  ADR_tem % kfl_laplacian          =  0 
  ADR_tem % kfl_nonlinear_sgs      =  kfl_sgsno_tem
  ADR_tem % kfl_time_sgs           =  kfl_sgsti_tem
  ADR_tem % kfl_time_bubble        =  kfl_tibub_tem
  ADR_tem % kfl_time_scheme        =  kfl_tisch_tem
  ADR_tem % kfl_time_order         =  kfl_tiacc_tem
  ADR_tem % kfl_manufactured       =  kfl_exacs_tem
  ADR_tem % kfl_length             =  kfl_ellen_tem
  ADR_tem % kfl_discretization     =  kfl_discr_tem
  if( kfl_sgsac_tem /= 1 ) then
     ADR_tem % kfl_first_order_sgs = 0
  else
     ADR_tem % kfl_first_order_sgs = 1
  end if
  ADR_tem % number_euler_steps     =  neule_tem

  ADR_tem % lun_output4            =  int(momod(modul) % lun_outpu,4)
  ADR_tem % bemol                  =  bemol_tem
  ADR_tem % tau_parameters(1:3)    =  staco_tem(1:3)
  ADR_tem % shock                  =  shock_tem  
  call ADR_check_and_compute_data(ADR_tem)
  call ADR_allocate_projections_bubble_sgs(ADR_tem)

  tesgs => ADR_tem % sgs

end subroutine tem_memall 
      
