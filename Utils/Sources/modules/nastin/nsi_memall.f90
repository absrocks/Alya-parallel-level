!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_memall.f90
!> @date    06/06/1966
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory
!> @details Allocate memory
!> @}
!-----------------------------------------------------------------------
subroutine nsi_memall()
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_nastin
  use mod_memory
  use def_kermod,              only : kfl_adj_prob,kfl_ndvars_opt,kfl_dvar_type,kfl_wlaav_ker, velav_ker,kfl_waexl_ker,avupo_ker
  use def_kermod,              only : kfl_twola_ker,kfl_noslw_ker
  use mod_communications,      only : PAR_MIN
  use mod_communications,      only : PAR_MAX
  use mod_parall,              only : PAR_GLOBAL_TO_LOCAL_NODE
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  implicit none
  integer(ip) :: ncsgs,ielem,pgaus,pelty,ipoin,pnode
  integer(ip) :: max_code,iboun,ivari

  !----------------------------------------------------------------------
  !
  ! Solver
  !
  !----------------------------------------------------------------------
  !
  ! Memory
  !
  solve_sol => solve(1:7)
  call soldef(4_ip)
  !
  ! Boundary conditions
  !
  solve(1) % bvess     => bvess_nsi(:,:,1)   ! Momentum
  solve(1) % kfl_fixno => kfl_fixno_nsi

  if( solve(1) % kfl_iffix /= 0 ) then       ! Velocity correction
     solve(4) % kfl_iffix =  2
     solve(4) % kfl_fixno => kfl_fixno_nsi
     solve(4) % bvess     => null()
  end if

  solve(2) % kfl_fixno => kfl_fixpr_nsi      ! Pressure schur complement
  solve(2) % bvess     => null()
  !
  ! Navier-Stokes is block-wise:
  ! First block:  momentum equations
  ! Second block: continuity equation
  !
  solve(1) % block_array(1) % bvess     => bvess_nsi(:,:,1)
  solve(1) % block_array(1) % kfl_fixno => kfl_fixno_nsi
  solve(1) % block_array(1) % bvnat     => solve(1) % bvnat

  solve(1) % block_array(2) % bvess     => bpess_nsi(:,:,1)
  solve(1) % block_array(2) % kfl_fixno => kfl_fixpp_nsi
  solve(1) % block_array(2) % bvnat     => solve(2) % bvnat
  !
  ! Tell the solver that although we have defined
  ! 2 solvers (momentum and continuity), the
  ! matrix used for these two solvers has a block
  ! structure, owned by the first solver!
  !
  if( NSI_MONOLITHIC ) then
     solve(1) % num_blocks          = 1
     solve(1) % block_dimensions(1) = ndime+1
  else if( NSI_SCHUR_COMPLEMENT ) then
     solve(1) % num_blocks          = 2
     solve(1) % block_dimensions(1) = ndime
     solve(1) % block_dimensions(2) = 1
     solve(1) % block_num           = 1
     solve(2) % block_num           = 2
  else if( NSI_FRACTIONAL_STEP ) then
     solve(1) % num_blocks          = 2
     solve(1) % block_dimensions(1) = ndime
     solve(1) % block_dimensions(2) = 1
     solve(1) % block_num           = 1
     solve(2) % block_num           = 2
  end if
  !
  ! Divergence free correction
  !
  if( kfl_divcorrec_nsi /= 0 ) then
     solve(6) % kfl_iffix =  2
     solve(6) % kfl_fixno => kfl_fixno_div_nsi
     solve(6) % bvess     => null()
  end if
  !
  ! Mass correction  - actually I guess that this lines would not be necesary if the dirichlet bcs are not imposed by the solver
  !
  if( kfl_corre_nsi == 3 ) then
     solve(4) % kfl_iffix =  0               ! not done by the solver
     solve(4) % bvess     => null()
     solve(4) % num_blocks          = 1
     solve(4) % block_dimensions(1) = ndime
  end if
  !
  ! Schur complement system: modify matrix and RHIS size
  ! Velocity and pressure are consecutively in UNKNO even
  ! if split scheme is used
  !
  if( INOTMASTER ) then

     if( NSI_SCHUR_COMPLEMENT ) then
        !
        ! Schur complement solver
        !
        nmauu_nsi =  solve(1) % nzmat
        nmaup_nsi =  solve(1) % nzmat/ndime
        nmapu_nsi =  solve(1) % nzmat/ndime
        nmapp_nsi =  solve(2) % nzmat
        poauu_nsi =  1
        poaup_nsi =  poauu_nsi + nmauu_nsi
        poapu_nsi =  poaup_nsi + nmaup_nsi
        poapp_nsi =  poapu_nsi + nmapu_nsi
        nzmat     =  max(nzmat,nmauu_nsi+nmaup_nsi+nmapu_nsi+nmapp_nsi)
        nzrhs     =  max(nzrhs,(ndime+1_ip)*npoin)

     else if( NSI_FRACTIONAL_STEP ) then
        !
        ! Fractional step: Auu and App are not assembles
        !
        nmauu_nsi =  solve(1) % nzmat
        nmaup_nsi =  solve(1) % nzmat/ndime
        nmapu_nsi =  solve(1) % nzmat/ndime
        nmapp_nsi =  solve(2) % nzmat
        poauu_nsi =  1
        poaup_nsi =  poauu_nsi + nmauu_nsi
        poapu_nsi =  poaup_nsi + nmaup_nsi
        poapp_nsi =  poapu_nsi + nmapu_nsi
        nzmat     =  max(nzmat,nmauu_nsi+nmaup_nsi+nmapu_nsi+nmapp_nsi)
        nzrhs     =  max(nzrhs,(ndime+1_ip)*npoin)
        !nmauu_nsi =  0
        !nmaup_nsi =  solve(1) % nzmat/ndime
        !nmapu_nsi =  solve(1) % nzmat/ndime
        !nmapp_nsi =  0
        !poauu_nsi =  1
        !poaup_nsi =  poauu_nsi + nmauu_nsi
        !poapu_nsi =  poaup_nsi + nmaup_nsi
        !poapp_nsi =  poapu_nsi + nmapu_nsi
        !nzmat     =  max(nzmat,nmauu_nsi+nmaup_nsi+nmapu_nsi+nmapp_nsi)
        !nzrhs     =  max(nzrhs,(ndime+1_ip)*npoin)
     end if

  else

     nmauu_nsi =  1
     nmaup_nsi =  1
     nmapu_nsi =  1
     nmapp_nsi =  1
     poauu_nsi =  1
     poaup_nsi =  1
     poapu_nsi =  1
     poapp_nsi =  1

  end if

  !----------------------------------------------------------------------
  !
  ! Arrays
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then
     !
     ! VELOC, PRESS: Allocate memory for velocity and pressure
     ! VELOC(:,:,1) = u^{n,i}
     ! VELOC(:,:,2) = u^{n,i-1}
     ! VELOC(:,:,3) = u^{n-1}
     ! VELOC(:,:,4) = u^{n-2}
     ! VELOC(:,:,5) = u^{n-3}, etc.
     !
     call memory_alloca(mem_modul(1:2,modul),'VELOC','nsi_memall',veloc,ndime,npoin,ncomp_nsi)
     call memory_alloca(mem_modul(1:2,modul),'PRESS','nsi_memall',press,npoin,ncomp_nsi)
     !
     ! DENSI: Compressible regime
     !
     if( kfl_regim_nsi == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'DENSI','nsi_memall',densi,npoin,1_ip)
     else if( kfl_regim_nsi == 2 ) then
        call memory_alloca(mem_modul(1:2,modul),'DENSI','nsi_memall',densi,npoin,ncomp_nsi)
     end if
     !
     ! Projetcions of dt/rho and tau
     !
     if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
        call memory_alloca(mem_modul(1:2,modul),'DT_RHO_NSI','nsi_memall',dt_rho_nsi,npoin)
        call memory_alloca(mem_modul(1:2,modul),'MU_RHO_NSI','nsi_memall',mu_rho_nsi,npoin)
        call memory_alloca(mem_modul(1:2,modul),'TAU_NSI'   ,'nsi_memall',tau_nsi   ,npoin)
     end if
     !
     ! Bubble
     !
     if( kfl_bubbl_nsi /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'BUBBLE_NSI'    ,'nsi_memall',bubble_nsi    ,nelem)
        call memory_alloca(mem_modul(1:2,modul),'BUBBLE_AQQ_NSI','nsi_memall',bubble_aqq_nsi,nelem)
        call memory_alloca(mem_modul(1:2,modul),'BUBBLE_AQU_NSI','nsi_memall',bubble_aqu_nsi,mnode*ndime,nelem)
        call memory_alloca(mem_modul(1:2,modul),'BUBBLE_AQP_NSI','nsi_memall',bubble_aqp_nsi,mnode,nelem)
        call memory_alloca(mem_modul(1:2,modul),'BUBBLE_BQ_NSI' ,'nsi_memall',bubble_bq_nsi ,nelem)
     end if
     !
     ! Immersed boundary method
     !
     if(      kfl_immer_nsi == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'LAGRA_NSI','nsi_memall',lagra_nsi,ndime,npoin,1_ip)
     else if( kfl_immer_nsi == 2 ) then
        call memory_alloca(mem_modul(1:2,modul),'TAUIB_NSI','nsi_memall',tauib_nsi,ndime,ndime,npoin)
     end if
     !
     ! UNK2N_NSI, Nastin second variable: PRESS or DENSI
     !
     if( kfl_regim_nsi == 2 ) then
        unk2n_nsi => densi
     else
        unk2n_nsi => press
     end if
     !
     ! VEOLD_NSI: Allocate memory for old (last iteration) velocity
     !
     ivari = 4
     call posdef(25_ip,ivari)
     if( ivari > 0 ) then
        kfl_resid_nsi=1
        call memory_alloca(mem_modul(1:2,modul),'VEOLD','nsi_memall',veold_nsi,ndime,npoin)
     else
        kfl_resid_nsi=0
     end if
     !
     ! VEPRO_NSI, PRPRO_NSI, GRPRO_NSI: Projections for orthogonal SGS
     !
     if( kfl_stabi_nsi > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'VEPRO_NSI','nsi_memall',vepro_nsi,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'PRPRO_NSI','nsi_memall',prpro_nsi,npoin)
        if( kfl_stabi_nsi == 2 ) then
           call memory_alloca(mem_modul(1:2,modul),'GRPRO_NSI','nsi_memall',grpro_nsi,ndime,npoin)
        end if
     end if
     !
     ! VESGS: Subgrid scale velocity
     !
     if( kfl_sgsco_nsi == 1 .or. kfl_sgsti_nsi == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'VESGS','nsi_memall',vesgs,nelem)
        ncsgs=min(2_ip,2_ip*kfl_sgsti_nsi+kfl_sgsco_nsi)
        do ielem=1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'VESGS(IELEM)','nsi_memall',vesgs(ielem)%a,ndime,pgaus,ncsgs)
        end do
     end if
     !
     ! DUNKN_NSI, DUNKP_NSI: Delta velocity and pressure  for Aitken relaxation strategy
     !
     if(kfl_relax_nsi==2) then
        call memory_alloca(mem_modul(1:2,modul),'DUNKN_NSI','nsi_memall',dunkn_nsi,ndime*npoin)
     end if
     if(kfl_relap_nsi==2) then
        call memory_alloca(mem_modul(1:2,modul),'DUNKN_NSI','nsi_memall',dunkn_nsi,npoin)
     end if
     !
     ! AVVEL_NSI: average velocity (postprocess)
     !
     if(postp(1) % npp_stepi(21)/=0.or.maxval(postp(1) % pos_times(1:nvart,21))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVVEL_NSI','nsi_memall',avvel_nsi,ndime,npoin)
     end if
     !
     ! AVPRE_NSI: average pressure (postprocess)
     !
     if(postp(1) % npp_stepi(28)/=0.or.maxval(postp(1) % pos_times(1:nvart,28))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVPRE_NSI','nsi_memall',avpre_nsi,npoin)
     end if
     !
     ! AVVE2_NSI: average velocity**2 (postprocess)
     !
     if(postp(1) % npp_stepi(57)/=0.or.maxval(postp(1) % pos_times(1:nvart,57))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVVE2_NSI','nsi_memall',avve2_nsi,ndime,npoin)
     end if
     !
     ! AVVXY_NSI: average vx*vy (postprocess)
     !
     if(postp(1) % npp_stepi(58)/=0.or.maxval(postp(1) % pos_times(1:nvart,58))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVVXY_NSI','nsi_memall',avvxy_nsi,ndime,npoin)
     end if
     !
     ! AVPRE_NSI: average pressure (postprocess)
     !
     if(postp(1) % npp_stepi(59)/=0.or.maxval(postp(1) % pos_times(1:nvart,59))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVPR2_NSI','nsi_memall',avpr2_nsi,npoin)
     end if
     !
     ! AVTAN_NSI: average TANGE (postprocess)
     !
     if(postp(1) % npp_stepi(46)/=0.or.maxval(postp(1) % pos_times(1:nvart,46))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVTAN_NSI','nsi_memall',avtan_nsi,ndime,npoin)
     end if
     !
     ! RESCH_NSI: Schur complement residual
     !
     if(postp(1) % npp_stepi(22)/=0.or.maxval(postp(1) % pos_times(1:nvart,22))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'RESCH_NSI','nsi_memall',resch_nsi,npoin)
     end if
     !
     ! REMOM_NSI: Schur complement residual
     !
     if(postp(1) % npp_stepi(62)/=0.or.maxval(postp(1) % pos_times(1:nvart,62))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'REMOM_NSI','nsi_memall',remom_nsi,ndime,npoin)
     end if
     !
     ! TURMU_NSI: LES turbulent viscosity (postprocess)
     !
     if(postp(1) % npp_stepi(74)/=0.or.maxval(postp(1) % pos_times(1:nvart,74))>zensi &
          .or. postp(1) % npp_stepi(75)/=0.or.maxval(postp(1) % pos_times(1:nvart,75))>zensi &
          .or. postp(1) % npp_stepi(76)/=0.or.maxval(postp(1) % pos_times(1:nvart,76))>zensi &
          .or. postp(1) % npp_stepi(77)/=0.or.maxval(postp(1) % pos_times(1:nvart,77))>zensi &
          .or. postp(1) % npp_stepi(78)/=0.or.maxval(postp(1) % pos_times(1:nvart,78))>zensi ) then
        call memory_alloca(mem_modul(1:2,modul),'TURMU_NSI','nsi_memall',turmu_nsi,nelem)
        do ielem=1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'TURMU_NSI','nsi_memall',turmu_nsi(ielem)%a,1_ip,pgaus,1_ip)
        end do
     end if
     !
     ! AVMUT_NSI: Average turbulent viscosity (postprocess)
     !
     if(postp(1) % npp_stepi(75)/=0.or.maxval(postp(1) % pos_times(1:nvart,75))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVMUT_NSI','nsi_memall',avmut_nsi,npoin)
     end if
     !
     ! AVSTX_NSI: Average stress mu_t * grad(u)
     !
     if(postp(1) % npp_stepi(76)/=0.or.maxval(postp(1) % pos_times(1:nvart,76))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVSTX_NSI','nsi_memall',avstx_nsi,ndime,npoin)
     end if
     !
     ! AVSTY_NSI: Average stress mu_t * grad(v)
     !
     if(postp(1) % npp_stepi(77)/=0.or.maxval(postp(1) % pos_times(1:nvart,77))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVSTY_NSI','nsi_memall',avsty_nsi,ndime,npoin)
     end if
     !
     ! AVSTZ_NSI: Average stress mu_t * grad(w)
     !
     if(postp(1) % npp_stepi(78)/=0.or.maxval(postp(1) % pos_times(1:nvart,78))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'AVSTZ_NSI','nsi_memall',avstz_nsi,ndime,npoin)
     end if
     !
     ! ENVEL_NSI: ensemble velocity (postprocess)
     !
     if(postp(1) % npp_stepi(86)/=0.or.maxval(postp(1) % pos_times(1:nvart,86))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENVEL_NSI','nsi_memall',envel_nsi,ndime,npoin)
     end if
     !
     ! ENPRE_NSI: ensemble pressure (postprocess)
     !
     if(postp(1) % npp_stepi(87)/=0.or.maxval(postp(1) % pos_times(1:nvart,87))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENPRE_NSI','nsi_memall',enpre_nsi,npoin)
     end if
     !
     ! ENVE2_NSI: ensemble velocity**2 (postprocess)
     !
     if(postp(1) % npp_stepi(88)/=0.or.maxval(postp(1) % pos_times(1:nvart,88))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENVE2_NSI','nsi_memall',enve2_nsi,ndime,npoin)
     end if
     !
     ! ENVXY_NSI: ensemble vx*vy (postprocess)
     !
     if(postp(1) % npp_stepi(89)/=0.or.maxval(postp(1) % pos_times(1:nvart,89))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENVXY_NSI','nsi_memall',envxy_nsi,ndime,npoin)
     end if
     !
     ! ENPRE_NSI: ensemble pressure (postprocess)
     !
     if(postp(1) % npp_stepi(90)/=0.or.maxval(postp(1) % pos_times(1:nvart,90))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENPR2_NSI','nsi_memall',enpr2_nsi,npoin)
     end if
     !
     ! ENTAN_NSI: ensemble TANGE (postprocess)
     !
     if(postp(1) % npp_stepi(91)/=0.or.maxval(postp(1) % pos_times(1:nvart,91))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENTAN_NSI','nsi_memall',entan_nsi,ndime,npoin)
     end if
     !
     ! ENMUT_NSI: Ensemble turbulent viscosity (postprocess)
     !
     if(postp(1) % npp_stepi(92)/=0.or.maxval(postp(1) % pos_times(1:nvart,92))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENMUT_NSI','nsi_memall',enmut_nsi,npoin)
     end if
     !
     ! ENSTX_NSI: Ensemble stress mu_t * grad(u)
     !
     if(postp(1) % npp_stepi(93)/=0.or.maxval(postp(1) % pos_times(1:nvart,93))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENSTX_NSI','nsi_memall',enstx_nsi,ndime,npoin)
     end if
     !
     ! ENSTY_NSI: Ensemble stress mu_t * grad(v)
     !
     if(postp(1) % npp_stepi(94)/=0.or.maxval(postp(1) % pos_times(1:nvart,94))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENSTY_NSI','nsi_memall',ensty_nsi,ndime,npoin)
     end if
     !
     ! ENSTZ_NSI: Ensemble stress mu_t * grad(w)
     !
     if(postp(1) % npp_stepi(95)/=0.or.maxval(postp(1) % pos_times(1:nvart,95))>zensi) then
        call memory_alloca(mem_modul(1:2,modul),'ENSTZ_NSI','nsi_memall',enstz_nsi,ndime,npoin)
     end if
     !
     ! Hydrostatic level set: used to compute rho_{hyd} to add to NS' RHS
     !
     if( kfl_hydro_gravity_nsi /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'HYDRO_DENSITY_NSI','nsi_memall',hydro_density_nsi,nelem)
        do ielem=1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           call memory_alloca(mem_modul(1:2,modul),'HYDRO_DENSITY_NSI(IELEM)','nsi_memall',hydro_density_nsi(ielem)%a,pgaus)
        end do
     end if
     !
     ! Surface tension
     !
     if( kfl_surte_nsi /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'NORLE_NSI','nsi_memall',norle_nsi,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'CURLE_NSI','nsi_memall',curle_nsi,npoin)
     end if
     !
     ! Traction on boundary nodes for alternative AVTAN postprocess
     !
!     if( output_postprocess_check_variable_postprocess(103_ip) .or. kfl_noslw_ker /= 0 )  & ! for the moment I am allocating always else I ahve to put ifs in other subros so taht it does not acces non allocated array
     call memory_alloca(mem_modul(1:2,modul),'NOTRA_NSI','nsi_memall',notra_nsi,ndime,npoin)
     if( kfl_waexl_ker /= 0 .or. kfl_twola_ker /= 0 .or. kfl_noslw_ker /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'MASSB_NSI','nsi_memall',massb_nsi,npoin)
     end if
     !
     ! arrays for no slip wall law - later they can be used for other cases too
     ! I have intialized tehm to 0 here but ther better option might be possible
     !
     if( output_postprocess_check_variable_postprocess(101_ip) .or. kfl_noslw_ker /= 0 )  then
        call memory_alloca(mem_modul(1:2,modul),'VAFOR_NSI','nsi_memall',vafor_nsi,ndime,npoin)
        vafor_nsi = 0.0_rp
     end if
     if( output_postprocess_check_variable_postprocess(102_ip) .or. kfl_noslw_ker /= 0 )  then
        call memory_alloca(mem_modul(1:2,modul),'AVVAF_NSI','nsi_memall',avvaf_nsi,ndime,npoin)
        avvaf_nsi = 0.0_rp
     end if

     if( output_postprocess_check_variable_postprocess(104_ip) .or. kfl_noslw_ker /= 0 )  then
        call memory_alloca(mem_modul(1:2,modul),'AVNTR_NSI','nsi_memall',avntr_nsi,ndime,npoin)
        avntr_nsi = 0.0_rp
     end if
     if( output_postprocess_check_variable_postprocess(105_ip) .or. kfl_noslw_ker /= 0 )  then
        call memory_alloca(mem_modul(1:2,modul),'AVGTR_NSI','nsi_memall',avgtr_nsi,ndime,npoin)
        avgtr_nsi = 0.0_rp
     end if
     !
     ! Forces due to porous media
     !
!     if( ???? /= 0 ) then   ! decidir bien cuando hay que alocarlo
        call memory_alloca(mem_modul(1:2,modul),'BUPOR_NSI','nsi_memall',bupor_nsi,ndime,npoin)
!     end if
     !
     ! Traction on boundary nodes calculated from an auxiliary RANS simulation
     !
     if( kfl_twola_ker /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'BTRAC_NSI','nsi_memall',btrac_nsi,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'TRACR_NSI','nsi_memall',tracr_nsi,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'TLUAV_NSI','nsi_memall',tluav_nsi,ndime,npoin)
     end if
     !
     ! Time-averaged velocity for wall law
     !
     if(kfl_wlaav_ker/=0) then
        call memory_alloca(mem_modul(1:2,modul),'VELAV_KER','nsi_memall',velav_ker,ndime,mgaub,max(1_ip,nboun))
!       if(kfl_waexl_ker==0) then
           call memory_alloca(mem_modul(1:2,modul),'AVUPO_KER','nsi_memall',avupo_ker,ndime,max(1_ip,npoin))
!       end if
     end if
     !
     ! Laplacian matrix
     !
     if( kfl_predi_nsi /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'LAPLA_NSI','nsi_memall',lapla_nsi,solve(2) % nzmat)
     end if
     !
     ! Consistent mass matrix
     !
     if(kfl_corre_nsi==3) then
        call memory_alloca(mem_modul(1:2,modul),'CMAMA_NSI','nsi_memall',cmama_nsi,solve(4) % nzmat)
        write(*,*) 'CONTROL LUEGO QUITAR: solve(2) % nzmat, solve(4) % nzmat',solve(2) % nzmat, solve(4) % nzmat
        if(solve(2) % nzmat*ndime*ndime /= solve(4) % nzmat) call runend('NSI_MEMALL:solve(2) % nzmat*ndime*ndime /= solve(4) % nzmat - PERHAPS you forgot MASSCORRECTION solver') !QUITARLO una vez probado
     end if
     !
     ! Save linear matrix
     !
     if(kfl_savco_nsi==1) then
        call memory_alloca(mem_modul(1:2,modul),'AMATR_NSI','nsi_memall',amatr_nsi,solve(1) % nzmat)
     end if
     !
     ! Coupling with SOLIDZ
     !
     if( coupling('SOLIDZ','NASTIN') >= 1 .or. coupling('NASTIN','IMMBOU') >= 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'FORCF','nsi_memall',forcf,ndime,npoin)
     end if
     !
     ! Matrix copy for internal force
     !
     if( kfl_intfo_nsi > 0 ) then
        allocate(intfo_nsi(npoin))
        do ipoin = 1,npoin
           intfo_nsi(ipoin) % kfl_exist = 0
        end do
     end if

     if (kfl_bnods_nsi == 1) then
        allocate(iboun_nsi(npoin))
     end if
     !
     ! Optimization materials for adjoint solution
     !
     if( kfl_adj_prob == 1 ) then
        if( kfl_dvar_type == 5 ) kfl_ndvars_opt = ndime*npoin
        call memory_alloca(mem_modul(1:2,modul),'VELOC_FORW','nsi_memall',veloc_forw,ndime,npoin,ncomp_nsi)
        call memory_alloca(mem_modul(1:2,modul),'PRESS_FORW','nsi_memall',press_forw,npoin,ncomp_nsi)
        !       call memory_alloca(mem_modul(1:2,modul),'RESDIFF_NSI','nsi_memall',resdiff_nsi,kfl_ndvars_opt, nzrhs)
        call memory_alloca(mem_modul(1:2,modul),'DCOST_DX_NSI','nsi_memall',dcost_dx_nsi,ndime*npoin)
        ! Coupling with temper
        if(kfl_coupl(ID_TEMPER,ID_NASTIN) == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'RhsadjTem_nsi','nsi_memall',RhsadjTem_nsi,nelem)
           do ielem=1,nelem
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              call memory_alloca(mem_modul(1:2,modul),'RhsadjTem_nsi','nsi_memall',RhsadjTem_nsi(ielem)%a,pnode)
           end do
        end if
        ! Coupling with turbul
        nturb = 2_ip
        if( kfl_coupl(ID_TURBUL,ID_NASTIN) == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'RhsadjTur_nsi','nsi_memall',RhsadjTur_nsi,nelem)
           do ielem=1,nelem
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              call memory_alloca(mem_modul(1:2,modul),'RhsadjTur_nsi','nsi_memall',RhsadjTur_nsi(ielem)%a,nturb,pnode)
           end do
        end if
     end if !adjoint

  else
     !
     ! Master: allocate minimum memory
     !
     call memory_alloca(mem_modul(1:2,modul),'VELOC','nsi_memall',veloc,1_ip,1_ip,ncomp_nsi)
     call memory_alloca(mem_modul(1:2,modul),'PRESS','nsi_memall',press,1_ip,ncomp_nsi)
     call memory_alloca_min(tluav_nsi)
     call memory_alloca_min(btrac_nsi)
     call memory_alloca_min(tracr_nsi)
     if( kfl_regim_nsi == 2 ) then
        unk2n_nsi => densi
     else
        unk2n_nsi => press
     end if
     if( kfl_stabi_nsi > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'VEPRO_NSI','nsi_memall',vepro_nsi,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'PRPRO_NSI','nsi_memall',prpro_nsi,1_ip)
        if( kfl_stabi_nsi == 2 ) then
           call memory_alloca(mem_modul(1:2,modul),'GRPRO_NSI','nsi_memall',grpro_nsi,1_ip,1_ip)
        end if
     end if
     if(kfl_relax_nsi==2) then
        call memory_alloca(mem_modul(1:2,modul),'DUNKN_NSI','nsi_memall',dunkn_nsi,1_ip)
     end if
     if(kfl_relap_nsi==2) then
        call memory_alloca(mem_modul(1:2,modul),'DUNKP_NSI','nsi_memall',dunkp_nsi,1_ip)
     end if
     if(kfl_regim_nsi==1.or.kfl_regim_nsi==2) then
        call memory_alloca(mem_modul(1:2,modul),'DENSI','nsi_memall',densi,1_ip,1_ip)
     end if
     if(kfl_predi_nsi/=0) then
        call memory_alloca(mem_modul(1:2,modul),'LAPLA_NSI','nsi_memall',lapla_nsi,1_ip)
     end if
     if(kfl_corre_nsi==3) then
        call memory_alloca(mem_modul(1:2,modul),'CMAMA_NSI','nsi_memall',cmama_nsi,1_ip)
     end if
     if( kfl_hydro_gravity_nsi /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'HYDRO_DENSITY_NSI','nsi_memall',hydro_density_nsi,1_ip)
     end if
     if( NSI_FRACTIONAL_STEP ) then
        call memory_alloca(mem_modul(1:2,modul),'DT_RHO_NSI','nsi_memall',dt_rho_nsi,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'MU_RHO_NSI','nsi_memall',mu_rho_nsi,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'TAU_NSI'   ,'nsi_memall',tau_nsi   ,1_ip)
     end if
     !
     ! Optimization materials for adjoint solution
     !
     if( kfl_adj_prob == 1 ) then
        if( kfl_dvar_type == 5 ) kfl_ndvars_opt = ndime*npoin
        call memory_alloca(mem_modul(1:2,modul),'DCOST_DX_NSI','nsi_memall',dcost_dx_nsi,1)
     end if !adjoint
     !
     ! VESGS: Subgrid scale velocity (needed for output of residuals)
     !
     if( kfl_sgsco_nsi == 1 .or. kfl_sgsti_nsi == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'VESGS','nsi_memall',vesgs,1_ip)
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! For condition of 20 type, need to know how much outflows we have
  !
  !----------------------------------------------------------------------

  if( kfl_exist_fib20_nsi /= 0 ) then
     max_code = -huge(0_ip)
     if (INOTMASTER) then
        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) == 20 ) max_code = max(max_code,int(bvnat_nsi(4,iboun,1),ip))
        end do
     end if

     call PAR_MAX(max_code)
     if( max_code < 1 ) call runend('NSI_MEMALL: WRONG DEFINITION OF STABLE OUTFLOW CONDITION')
     call memory_alloca(mem_modul(1:2,modul),'OUTFLOW_MASS','nsi_memall',outflow_mass,max_code)
  end if

end subroutine nsi_memall
