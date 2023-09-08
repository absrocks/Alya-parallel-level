subroutine ker_parall(order)
  !-----------------------------------------------------------------------
  !****f* Parall/ker_parall
  ! NAME
  !    ker_parall
  ! DESCRIPTION
  !    This routine exchange data
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_inpout
  use def_solver
  use mod_memory
  use mod_opebcs
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  use mod_communications,     only : PAR_BROADCAST
  use mod_communications,     only : PAR_EXCHANGE
  use mod_ker_regularization, only : kfl_regularization, kfl_second, reg_type
  use mod_ker_subdomain,      only : ker_subdomain_parall
  use mod_ker_tendencies,     only : kfl_tendencies_ker
  implicit none
  integer(ip), intent(in)   :: order
  integer(ip)               :: ii,ji,dummi,imate,ipala,ki,ifunc
  integer(ip)               :: idime,kdime,nexpr
  
  if( ISEQUEN ) return

  strre = 'ker_paral'
  strin = 'ker_paral'
  strch = 'ker_paral'
  
  nullify(parin)
  nullify(parre)
  nullify(parlo)
  nullify(parhh)

  select case ( order )

  case( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Broadcast data read in *.ker.dat file
     !
     !-------------------------------------------------------------------
     
     do parii = 1,2
        
        npari = 0
        nparr = 0
        nparl = 0
        nparc = 0
        nparh = 0
        !
        ! Physical problem
        !
        call iexcha(kfl_vefun)
        call iexcha(kfl_rough)
        call iexcha(kfl_canhe)
        call iexcha(kfl_heiov)
        call iexcha(kfl_canla)
        call iexcha(kfl_delta)
        call iexcha(kfl_ustar)
        call iexcha(kfl_walld)
        call iexcha(kfl_walln)
        do ji = 1,2
           call iexcha(kfl_walld_field(ji))
        end do
        call iexcha(kfl_suppo)
        call iexcha(kfl_extro)
        call iexcha(kfl_prope)
        call iexcha(kfl_kemod_ker)
        call iexcha(kfl_noslw_ker)
        call iexcha(kfl_waexl_ker)
        call iexcha(kfl_waexl_imp_ker)
        call iexcha(kfl_twola_ker)
        do ji = 1,9
           call iexcha(krestr_2_codno(ji))
        end do
        call iexcha(nrestr_2_codno)
        call iexcha(kfl_wlaav_ker)
        call iexcha(kfl_aveme_ker)
        call iexcha(kfl_algebra_operations)
        do ji = 1,mcodb
           call iexcha(kfl_boexc_ker(ji))
        end do
        call iexcha(kfl_conta)

        call rexcha(denme)
        call rexcha(visme)
        call rexcha(gasco)
        call rexcha(u_ref)
        call rexcha(h_ref)
        call rexcha(k_ref)
        call rexcha(usref)

        do ji = 1,3
           call rexcha(windg(ji))
        end do

        call rexcha(delta_dom)
        call rexcha(delmu_dom)
        call rexcha(rough_dom)
        call rexcha(grnor)
        call rexcha(gravi(1))
        call rexcha(gravi(2))
        call rexcha(gravi(3))
        call rexcha(thicl)
        call rexcha(cmu_st)
        call rexcha(dexlo_ker)
        call rexcha(tpeav_ker)
        call iexcha(kfl_tendencies_ker)
        call rexcha(tlape_ker)
        call rexcha(cpres_ker(1))
        call rexcha(rpres_ker(1))
        call rexcha(cpres_ker(2))
        call rexcha(rpres_ker(2))
        !
        ! Coupling between modules
        !
        do ji = 0,mmodu
           do ki = 0,mmodu
              call iexcha(kfl_coupl(ki,ji))
              call iexcha(kfl_cowhe(ki,ji))
           end do
        end do
        !
        ! Properties
        !
        call iexcha(densi_ker % kfl_exist)
        call iexcha(visco_ker % kfl_exist)
        call iexcha(poros_ker % kfl_exist)
        call iexcha(condu_ker % kfl_exist)
        call iexcha(sphea_ker % kfl_exist)
        call iexcha(dummy_ker % kfl_exist)
        call iexcha(turmu_ker % kfl_exist)
        call iexcha(absor_ker % kfl_exist)
        call iexcha(scatt_ker % kfl_exist)

        do imate = 1,nmate
           do ipala = 1,mlapa_ker
              call rexcha(densi_ker % rlaws(ipala,imate))
              call rexcha(visco_ker % rlaws(ipala,imate))
              call rexcha(poros_ker % rlaws(ipala,imate))
              call rexcha(condu_ker % rlaws(ipala,imate))
              call rexcha(sphea_ker % rlaws(ipala,imate))
              call rexcha(dummy_ker % rlaws(ipala,imate))
              call rexcha(turmu_ker % rlaws(ipala,imate))
              call rexcha(absor_ker % rlaws(ipala,imate))
              call rexcha(scatt_ker % rlaws(ipala,imate))
           end do
        end do

        do imate = 1,nmate
           call rexcha(densi_ker % value_default(imate))
           call rexcha(visco_ker % value_default(imate))
           call rexcha(poros_ker % value_default(imate))
           call rexcha(condu_ker % value_default(imate))
           call rexcha(sphea_ker % value_default(imate))
           call rexcha(dummy_ker % value_default(imate))
           call rexcha(turmu_ker % value_default(imate))
           call rexcha(absor_ker % value_default(imate))
           call rexcha(scatt_ker % value_default(imate))
        end do

        do imate = 1,nmate
           call PAR_EXCHANGE(5_ip,densi_ker % wlaws(imate),parhh,nparh,parii) 
           call PAR_EXCHANGE(5_ip,visco_ker % wlaws(imate),parhh,nparh,parii) 
           call PAR_EXCHANGE(5_ip,poros_ker % wlaws(imate),parhh,nparh,parii) 
           call PAR_EXCHANGE(5_ip,condu_ker % wlaws(imate),parhh,nparh,parii) 
           call PAR_EXCHANGE(5_ip,sphea_ker % wlaws(imate),parhh,nparh,parii) 
           call PAR_EXCHANGE(5_ip,dummy_ker % wlaws(imate),parhh,nparh,parii) 
           call PAR_EXCHANGE(5_ip,turmu_ker % wlaws(imate),parhh,nparh,parii) 
           call PAR_EXCHANGE(5_ip,absor_ker % wlaws(imate),parhh,nparh,parii) 
           call PAR_EXCHANGE(5_ip,scatt_ker % wlaws(imate),parhh,nparh,parii) 
        end do

        do imate = 1,nmate
           do ii = 1, 2
              call iexcha(densi_ker % update(ii,imate))
              call iexcha(visco_ker % update(ii,imate))
              call iexcha(poros_ker % update(ii,imate))
              call iexcha(condu_ker % update(ii,imate))
              call iexcha(sphea_ker % update(ii,imate))
              call iexcha(dummy_ker % update(ii,imate))
              call iexcha(turmu_ker % update(ii,imate))
              call iexcha(absor_ker % update(ii,imate))
              call iexcha(scatt_ker % update(ii,imate))
           end do
        end do
        call ker_subdomain_parall()
        !
        ! Numerical problem
        !
        call iexcha(kfl_renumbering_npoin)
        call iexcha(kfl_renumbering_nelem)
        call iexcha(nsfc_renumbering_npoin)
        call iexcha(ndivi)
        call iexcha(multiply_with_curvature)
        call iexcha(kfl_edge_elements)
        call iexcha(kfl_rotation_axe)
        call iexcha(kfl_graph)
        call iexcha(kfl_elm_graph)
        call iexcha(kfl_savda)
        call iexcha(kfl_data_base_array)
        call iexcha(kfl_vector_size)
        call iexcha(kfl_lface)
        call iexcha(kfl_fv_data)
        call iexcha(kfl_grpro)
        call iexcha(kfl_fixsm)
        call iexcha(kfl_matrix_grad)
        call iexcha(kfl_conbc_ker)
        call iexcha(kfl_element_bin)
        call iexcha(kfl_elses)
        do ji=1,nelse
           call iexcha(ielse(ji))
        end do
        do ji=1,nelse
           call rexcha(relse(ji))
        end do
        call rexcha(rotation_angle)

        call iexcha(kfl_cutel)
        call iexcha(kfl_hangi)
        call iexcha(kfl_lapla)
        call iexcha(kfl_defor)
        call iexcha(kfL_coo)
        call iexcha(kfL_ell)
        call iexcha(kfl_full_rows)
        call iexcha(kfl_element_to_csr)
        call iexcha(kfl_direct_solver)
        call iexcha(npoin_mm)
        call iexcha(nboun_mm)
        call iexcha(mnodb_mm)

        call iexcha(number_space_time_function)                ! Space/Time functions
        do ifunc = 1,max_space_time_function
           call iexcha(space_time_function(ifunc) % ndime)
           call iexcha(space_time_function(ifunc) % nexpr)
           call cexcha(5_ip,space_time_function(ifunc) % name)
        end do
        call iexcha(number_time_function)                      ! Time functions
        do ifunc = 1,max_time_function
           call iexcha(time_function(ifunc) % npara)
           call iexcha(time_function(ifunc) % kfl_type)
           call cexcha(5_ip,time_function(ifunc) % name)
        end do

        call iexcha(deformation_steps)
        call iexcha(deformation_strategy)
        call iexcha(kfl_duatss)
        call iexcha(fact_duatss)
        call iexcha(kfl_conma)
        call iexcha(kfl_conma_weighted)
        call iexcha(kfl_logva)
        call iexcha(reg_type)
        !
        ! Output
        !
        call iexcha(nwitn)
        call iexcha(mwitn)
        call iexcha(kfl_posdo)
        call iexcha(kfl_posdi)
        do ii = 1,size(kfl_oumes)
           call iexcha(kfl_oumes(ii))
        end do
        call iexcha(kfl_oustl)
        call iexcha(kfl_quali)
        call iexcha(npp_stepo)
        call iexcha(nfilt)
        call iexcha(kfl_abovx)
        do ii = 1,3
           call iexcha(resvx(ii))
           call rexcha(travx(ii))
        end do
        call iexcha(kfl_vortx)
        call iexcha(kfl_vortx_thres)
        call iexcha(kfl_detection)
        call iexcha(kfl_pixel)
        call iexcha(kfl_livee)
        call iexcha(plane_pixel)
        call iexcha(variable_pixel)
        call iexcha(number_pixel(1))
        call iexcha(number_pixel(2))
        call iexcha(nsteps_ensemble)

        do ii = 1,2
           do ji = 1,3
              call rexcha(bobvx(ii,ji))
           end do
        end do
        do ii = 1,3
           call rexcha(travx(ii))
        end do
        call rexcha(thr_veloc)
        call rexcha(thr_qvort)
        call rexcha(detection_length)
        call rexcha(detection_velocity)
        call rexcha(coord_plane_pixel)

        solve_sol => solve(1:)
        call soldef(1_ip)        ! Solvers
        call fildef(1_ip)        ! Filters
        call posdef(1_ip,dummi)  ! Postprocess
	!
	! adjoint and optimization
	!
        call iexcha(kfl_cost_type)
        call iexcha(kfl_dvar_type)
        call iexcha(kfl_adj_prob)
        call iexcha(kfl_cos_opt)
        call iexcha(kfl_ndvars_opt)
        call iexcha(kfl_nwall)
        !
        ! Others
        ! ESTA LIENA ESTA MUY MAL!!!!!
        !call cexcha(len(events_directory),events_directory) ! Events directory name
        if( parii == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'PARIN','ker_parall',parin,npari,'DO_NOT_INITIALIZE')
           call memory_alloca(mem_modul(1:2,modul),'PARRE','ker_parall',parre,nparr,'DO_NOT_INITIALIZE')
           call memory_alloca(mem_modul(1:2,modul),'PARLO','ker_parall',parlo,nparl,'DO_NOT_INITIALIZE')
           call memory_alloca(mem_modul(1:2,modul),'PARHH','ker_parall',parhh,nparh,'DO_NOT_INITIALIZE')
        end if
        if( ( parii == 1 .and. ISLAVE ) .or. ( parii == 2 .and. IMASTER ) ) then
           call PAR_BROADCAST(parin,      'IN MY CODE')
           call PAR_BROADCAST(parre,      'IN MY CODE')
           call PAR_BROADCAST(parlo,      'IN MY CODE')
           call PAR_BROADCAST(nparc,parch,'IN MY CODE')
           call PAR_BROADCAST(nparh,parhh,'IN MY CODE')
        end if
     end do

     call memory_deallo(mem_modul(1:2,modul),'PARIN','ker_parall',parin)
     call memory_deallo(mem_modul(1:2,modul),'PARRE','ker_parall',parre)
     call memory_deallo(mem_modul(1:2,modul),'PARLO','ker_parall',parlo)
     call memory_deallo(mem_modul(1:2,modul),'PARHH','ker_parall',parhh)
     !
     ! Broadcasting of logical data
     !
     call par_broadcast(kfl_regularization)
     call par_broadcast(kfl_second)

     !-------------------------------------------------------------------
     !
     ! Witness point coordinates
     !
     !-------------------------------------------------------------------

     if( nwitn /= 0 ) then
        if( ISLAVE ) call memose(7_ip)
        !call PAR_BROADCAST(nwitn*3_ip,cowit)
        call pararr('BCT',0_ip,nwitn*3,cowit)
     end if

     !-------------------------------------------------------------------
     !
     ! Boundary conditions
     !
     !-------------------------------------------------------------------

     call spnbcs(tncod_ker)
     call spgbcs(tgcod_ker)
     call spbbcs(tbcod_ker)

     !-------------------------------------------------------------------
     !
     ! Space/Time functions
     !
     !-------------------------------------------------------------------

     if( number_space_time_function > 0 ) then
        if( ISLAVE ) then
           do ifunc = 1,number_space_time_function
              igene = ifunc
              call ker_memory(7_ip)
           end do
        end if
        do parii = 1,2
           npari = 0
           nparr = 0
           nparc = 0
           do ifunc = 1,number_space_time_function
              nexpr = space_time_function(ifunc) % nexpr
              idime = space_time_function(ifunc) % ndime
              do kdime = 1,idime
                 call cexcha(nexpr,space_time_function(ifunc) % expression(kdime))
              end do
           end do
           if( parii == 1 ) then
              call memory_alloca(mem_modul(1:2,modul),'PARIN','ker_parall',parin,npari,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_modul(1:2,modul),'PARRE','ker_parall',parre,nparr,'DO_NOT_INITIALIZE')
              if( ISLAVE ) call Parall(2_ip)
           end if
        end do
        if( IMASTER ) call Parall(2_ip)
        call memory_deallo(mem_modul(1:2,modul),'PARRE','ker_parall',parre)
        call memory_deallo(mem_modul(1:2,modul),'PARIN','ker_parall',parin)
        if( IMASTER ) then
           do ifunc = 1,number_space_time_function
              igene = ifunc
              call ker_memory(-7_ip)
           end do
        end if
     end if

     !-------------------------------------------------------------------
     !
     ! Time functions
     !
     !-------------------------------------------------------------------

     if( number_time_function > 0 ) then
        if( ISLAVE ) then
           do ifunc = 1,number_time_function
              igene = ifunc
              call ker_memory(8_ip)
           end do
        end if
        do parii = 1,2
           npari = 0
           nparr = 0
           nparc = 0
           do ifunc = 1,number_time_function
              idime = time_function(ifunc) % npara
              do kdime = 1,idime
                 call rexcha(time_function(ifunc) % parameters(kdime))
              end do
           end do
           if( parii == 1 ) then
              call memory_alloca(mem_modul(1:2,modul),'PARIN','ker_parall',parin,npari,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_modul(1:2,modul),'PARRE','ker_parall',parre,nparr,'DO_NOT_INITIALIZE')
              if( ISLAVE ) call Parall(2_ip)
           end if
        end do

        if( IMASTER ) call Parall(2_ip)
        call memory_deallo(mem_modul(1:2,modul),'PARRE','ker_parall',parre)
        call memory_deallo(mem_modul(1:2,modul),'PARIN','ker_parall',parin)
        if( IMASTER ) then
           do ifunc = 1,number_time_function
              igene = ifunc
              call ker_memory(-8_ip)
           end do
        end if
     end if

  case( 2_ip )

     !-------------------------------------------------------------------
     !
     ! Support surface: master sends to slaves their corresponding parts
     !
     !-------------------------------------------------------------------

     strre = 'ker_paral'
     strin = 'ker_paral'
     strch = 'ker_paral'
     nullify(parin)
     nullify(parre)

     if( kfl_suppo == 1 ) then

        call parari('BCT',0_ip,1_ip,npoin_mm)
        call parari('BCT',0_ip,1_ip,nboun_mm)
        if( ISLAVE )  call ker_memory(5_ip)
        call pararr('BCT',0_ip,npoin_mm*ndime,   coord_mm)
        call parari('BCT',0_ip,nboun_mm*mnodb_mm,lnodb_mm)
        call parari('BCT',0_ip,nboun_mm         ,ltypb_mm)
        if( IMASTER ) call ker_memory(-5_ip)

!!$        if( ISLAVE ) then
!!$           !
!!$           ! Slave receive their corresponding geometry. Eventhough they do not receive
!!$           ! any boundary, NPOIN_MM should be different from zero to enter the loops
!!$           ! involving the support geometry treatment
!!$           !
!!$           kfl_desti_par = 0
!!$           call parari('RCV',0_ip,1_ip,npoin_mm)
!!$           call parari('RCV',0_ip,1_ip,nboun_mm)
!!$           call ker_memory(5_ip)
!!$
!!$           if( npoin_mm > 0 ) then
!!$              call pararr('RCV',0_ip,npoin_mm*ndime,   coord_mm)
!!$              call parari('RCV',0_ip,nboun_mm*mnodb_mm,lnodb_mm)
!!$              call parari('RCV',0_ip,nboun_mm         ,ltypb_mm)
!!$           end if
!!$
!!$        else
!!$           !
!!$           ! Send the support geometry inside in the slave's bounding box (and a bit more...)
!!$           !
!!$           call memory_alloca(mem_modul(1:2,modul),'XBOUN','ker_parall',xboun,2_ip,3_ip,nboun_mm)
!!$           call memory_alloca(mem_modul(1:2,modul),'GPOIN','ker_parall',gpoin,npoin_mm)
!!$           call memory_alloca(mem_modul(1:2,modul),'GBOUN','ker_parall',gboun,nboun_mm)
!!$
!!$           do iboun = 1,nboun_mm
!!$              xboun(1,1:3,iboun) =  1.0e12_rp
!!$              xboun(2,1:3,iboun) = -1.0e12_rp
!!$              do inodb = 1,mnodb_mm
!!$                 ipoin = lnodb_mm(inodb,iboun)
!!$                 do idime = 1,ndime
!!$                    if( coord_mm(idime,ipoin) < xboun(1,idime,iboun) ) xboun(1,idime,iboun) = coord_mm(idime,ipoin)
!!$                    if( coord_mm(idime,ipoin) > xboun(2,idime,iboun) ) xboun(2,idime,iboun) = coord_mm(idime,ipoin)
!!$                 end do
!!$              end do
!!$           end do
!!$
!!$           tomin = 0.9_rp
!!$           tomax = 1.1_rp
!!$
!!$           do ipart = 1,npart
!!$
!!$              kfl_desti_par = ipart
!!$              nboun_loc     = 0
!!$              npoin_loc     = 0
!!$              do iboun = 1,nboun_mm
!!$                 nboun_loc    = nboun_loc + 1
!!$                 gboun(iboun) = nboun_loc
!!$                 do inodb = 1,nnode(abs(ltypb_mm(iboun)))
!!$                    ipoin = lnodb_mm(inodb,iboun)
!!$                    if( gpoin(ipoin) == 0 ) then
!!$                       npoin_loc = npoin_loc + 1
!!$                       gpoin(ipoin) = npoin_loc
!!$                    end if
!!$                 end do
!!$              end do
!!$
!!$              if( npoin_loc > 0 ) then
!!$
!!$                 call memory_alloca(mem_modul(1:2,modul),'COORD_LOC','ker_parall',coord_loc,ndime,   npoin_loc,'DO_NOT_INITIALIZE')
!!$                 call memory_alloca(mem_modul(1:2,modul),'LNODB_LOC','ker_parall',lnodb_loc,mnodb_mm,nboun_loc,'DO_NOT_INITIALIZE')
!!$                 call memory_alloca(mem_modul(1:2,modul),'LTYPB_LOC','ker_parall',ltypb_loc,nboun_loc         ,'DO_NOT_INITIALIZE')
!!$
!!$                 do ipoin = 1,npoin_mm
!!$                    if( gpoin(ipoin) /= 0 ) then
!!$                       kpoin = gpoin(ipoin)
!!$                       do idime = 1,ndime
!!$                          coord_loc(idime,kpoin) = coord_mm(idime,ipoin)
!!$                       end do
!!$                    end if
!!$                 end do
!!$
!!$                 do iboun = 1,nboun_mm
!!$                    if( gboun(iboun) /= 0 ) then
!!$                       kboun = gboun(iboun)
!!$                       ltypb_loc(kboun) = ltypb_mm(iboun)
!!$                       do inodb = 1,mnodb_mm
!!$                          ipoin = lnodb_mm(inodb,iboun)
!!$                          lnodb_loc(inodb,kboun) = gpoin(ipoin)
!!$                       end do
!!$                       gboun(iboun) = 0
!!$                    end if
!!$                 end do
!!$                 do ipoin = 1,npoin_mm
!!$                    gpoin(ipoin) = 0
!!$                 end do
!!$
!!$              end if
!!$
!!$              call parari('SND',0_ip,1_ip,npoin_loc)
!!$              call parari('SND',0_ip,1_ip,nboun_loc)
!!$
!!$              if( npoin_loc > 0 ) then
!!$                 call pararr('SND',0_ip,npoin_loc*ndime,   coord_loc)
!!$                 call parari('SND',0_ip,nboun_loc*mnodb_mm,lnodb_loc)
!!$                 call parari('SND',0_ip,nboun_loc         ,ltypb_loc)
!!$                 call memory_deallo(mem_modul(1:2,modul),'LTYPB_LOC','ker_parall',ltypb_loc)
!!$                 call memory_deallo(mem_modul(1:2,modul),'LNODB_LOC','ker_parall',lnodb_loc)
!!$                 call memory_deallo(mem_modul(1:2,modul),'COORD_LOC','ker_parall',coord_loc)
!!$              end if
!!$
!!$           end do
!!$
!!$           call memory_deallo(mem_modul(1:2,modul),'GBOUN','ker_parall',gboun)
!!$           call memory_deallo(mem_modul(1:2,modul),'GPOIN','ker_parall',gpoin)
!!$           call memory_deallo(mem_modul(1:2,modul),'XBOUN','ker_parall',xboun)
!!$           !
!!$           ! Deallocate COORD_MM and LNODB_MM
!!$           !
!!$           call ker_memory(-5_ip)
!!$
!!$        end if
     end if

  end select

  npari = 0
  nparr = 0
  nparc = 0
  nullify(parin)
  nullify(parre)

end subroutine ker_parall

