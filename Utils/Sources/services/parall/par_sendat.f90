subroutine par_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Parall/par_sendat
  ! NAME
  !    par_sendat
  ! DESCRIPTION
  !    This routine exchange data
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_parall
  use def_inpout
  use def_solver
  use mod_memory
  use mod_parall
  use def_kermod
  use def_mpio
  use mod_memory, only : kfl_memor
  use mod_memory, only : kfl_varcount
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ii,ir,ic,ji,ki,dummi
  integer(ip)             :: ipass,kfl_ptask_old

  select case (order)

  case(1_ip)
     !
     ! First communication always performed with MPI and not files
     !
     if( IMASTER ) then
        if( kfl_ptask == 0 .and. nproc_par == 1 ) return
        kfl_ptask_old=kfl_ptask
     end if
     kfl_ptask=1
     call vocabu(-1_ip,0_ip,0_ip)
     !
     ! Exchange data read in Reapro (always using MPI)
     !
     if( IPARALL ) then

        strre='Reapro'
        strin='Reapro'
        strch='Reapro'
        do parii=1,2
           npari=0
           nparr=0
           nparc=0
           !
           ! Partition (just to check errors)
           !
           call iexcha(npart_par)
           !
           ! Read in rrudat
           !
           call iexcha(current_code)
           call iexcha(kfl_naked)
           call iexcha(kfl_custo)
           call iexcha(kfl_preli)
           call iexcha(kfl_memor)
           call iexcha(kfl_varcount)
           call iexcha(lun_livei)
           call iexcha(kfl_timin)
           call iexcha(kfl_lotme)
           call iexcha(kfl_freme)
           call iexcha(kfl_outpu)
           call iexcha(nprit)
           call iexcha(kfl_rstar)
           call iexcha(kfl_rsfil)
           call iexcha(kfl_timeline)
           call iexcha(kfl_commu)
           call iexcha(kfl_outfo)
           call iexcha(kfl_ptask_old)
           call iexcha(kfl_vtk)
           call iexcha(kfl_rread) ! Reread ini in READ_AND_RUN mode
#ifdef PARMETIS
           call iexcha(kfl_paral_parmes)
           call iexcha(kfl_paral_proc_node)
#endif
!           call iexcha(par_hybrid)
           call iexcha(par_omp_granularity)
           call iexcha(par_omp_coloring_alg)
           call iexcha(par_omp_nelem_chunk)
           call iexcha(par_omp_nboun_chunk)
           call iexcha(par_omp_npoin_chunk)
           call iexcha(par_omp_partition_alg)

           call iexcha(par_topo_num_nodes)
           call iexcha(par_topo_num_cores_per_node)

           call iexcha(kfl_matri_par)

           call iexcha(kfl_partition_par)
           call iexcha(kfl_parseq_par)
           call iexcha(kfl_interface_parti_par)
           call iexcha(kfl_weigh_par)
           call iexcha(kfl_parti_par)
           call iexcha(kfl_global_numbering_par)

           call iexcha(kfl_cores_per_gpu)
           call iexcha(kfl_streams_per_gpu)

           do ji = 1,size(weights_elements_par,KIND=ip)
              call rexcha(weights_elements_par(ji))
           end do
           do ji = 1,size(weights_materials_par,KIND=ip)
              call rexcha(weights_materials_par(ji))
           end do
           do ji = 1,3
              call iexcha(boxes_coarse_par(ji))
           end do
           do ji = 1,3
              call iexcha(boxes_fine_par(ji))
           end do
           do ji = 1,3
              call rexcha(vect_partition_par(ji))
           end do

           call iexcha(sfc_check)
           call iexcha(sfc_criteria)
           call iexcha(sfc_dim_bin_core)

           call rexcha(cpu_limit)
           nparc=nparc+132
           if(parii==2.and.IMASTER) parch(1:66)   = title(1:66)
           if(parii==2.and.ISLAVE)  title(1:66)   = parch(1:66)
           if(parii==2.and.IMASTER) parch(67:132) = namda(1:66)
           if(parii==2.and.ISLAVE)  namda(1:66)   = parch(67:132)
           call cexcha(int(len(method_redistribution_par),ip),method_redistribution_par)
           !
           ! Read in readat and 'mod'_reapro of modules 'mod'
           !
           do ji=1,mblok
              call iexcha(micou(ji))
           end do
           call iexcha(nblok)
           call iexcha(kfl_timco)
           call iexcha(kfl_reset)
           call rexcha(reset_factor)
           call iexcha(kfl_timei)
           call iexcha(kfl_timef)
           call iexcha(kfl_dtfun)
           call iexcha(mitim)
           call iexcha(mitsm)
           call iexcha(mitrf)
           call iexcha(kfl_algor_msh)
           call iexcha(kfl_error_msh)
           call iexcha(kfl_gover_msh)
           call iexcha(kfl_block_msh)
           call iexcha(kfl_outpu_par)
           call iexcha(kfl_postp_par)
           call iexcha(nzone_par)
           call iexcha(kfl_wwork)
           call iexcha(kfl_lumped)
           call rexcha(timei)
           call rexcha(timef)
           call rexcha(dtime)
           call rexcha(toler_msh)
           do ji=1,mmodu
              call iexcha(kfl_modul(ji))
           end do
           do ji=1,mmodu
              call iexcha(kfl_delay(ji))
           end do
           do ji=1,mmodu
              call iexcha(kfl_conve(ji))
           end do
           do ji=0,mmodu
              call iexcha(kfl_solve(ji))
           end do
           do ji=1,mmodu
              call iexcha(ndela(ji))
           end do
           do ji=1,mmodu
              do ki=1,mblok
                 call iexcha(lmord(ji,ki))
              end do
           end do
           do ji=0,mmodu
              call iexcha(lzone(ji))
           end do
           do ji=1,mserv
              call iexcha(kfl_servi(ji))
           end do
           call fildef(1_ip)
           call posdef(1_ip,dummi)
           !
           ! Parallel IO
           !
           call iexcha(mpio_flag_enabled)
           call iexcha(mpio_flag_geometry)
           call iexcha(mpio_flag_geometry_export)
           call iexcha(mpio_flag_geometry_read_post)
           call iexcha(mpio_flag_post)
           call iexcha(mpio_flag_post_light)
           call iexcha(mpio_flag_rst)
           call iexcha(mpio_flag_autoromio)
           call iexcha(mpio_flag_synchro)
           call iexcha(mpio_flag_collective)
           call iexcha(mpio_flag_communicator)
           call iexcha(mpio_flag_all_par)
           call iexcha(mpio_flag_post_merge)
           call iexcha(mpio_val_asyncbuffer)
           call iexcha(mpio_val_merge_block)
           call rexcha(mpio_val_hybrid_threshold)
           !
           ! Allocate memory for the first pass
           !
           if(parii==1) then
              call memory_alloca(mem_servi(1:2,servi),'PARIN','par_sendat',parin,npari)
              call memory_alloca(mem_servi(1:2,servi),'PARRE','par_sendat',parre,nparr)
              if(ISLAVE) call par_broadc()
           end if
        end do
        if(IMASTER) call par_broadc()

        call memory_deallo(mem_servi(1:2,servi),'PARIN','par_sendat',parin)
        call memory_deallo(mem_servi(1:2,servi),'PARRE','par_sendat',parre)

        kfl_ptask = kfl_ptask_old
        call vocabu(-1_ip,0_ip,0_ip)

     end if
     
  case(2_ip)

    if( IPARALL ) then

        strre='readim_reastr_reageo'
        strin='readim_reastr_reageo'
        strch='readim_reastr_reageo'
        do parii=1,2
           ipass = parii
           npari=0
           nparr=0
           nparc=0
           !
           ! Read in reageo
           !
           continue
           !
           ! Read in reaset
           !
           call iexcha(neset)
           call iexcha(nbset)
           call iexcha(nnset)     ! Re-computed in par_senset
           !
           ! Read in reabcs
           !
           call iexcha(kfl_icodn)
           call iexcha(kfl_icodb)
           !
           ! Read in reafie
           !
           continue

           if( ipass == 1 ) then
              call memory_alloca(mem_servi(1:2,servi),'PARIN','par_sendat',parin,npari)
              call memory_alloca(mem_servi(1:2,servi),'PARRE','par_sendat',parre,nparr)
              if( ISLAVE ) call par_receiv()
           end if
        end do

        if( IMASTER ) call par_sendin()

        call memory_deallo(mem_servi(1:2,servi),'PARIN','par_sendat',parin)
        call memory_deallo(mem_servi(1:2,servi),'PARRE','par_sendat',parre)

     end if

  case(5)

     if( kfl_paral >= 0 ) then

        strre='par_partit'
        strin='par_partit'
        strch='par_partit'
        do ipass=1,2
           ii=0
           ir=0
           ic=0
           !
           ! Calculated in partit
           !
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = gni
           if(ipass==2.and.ISLAVE)  gni         = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = gnb
           if(ipass==2.and.ISLAVE)  gnb         = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = ginde_par(3,kfl_desti_par)
           if(ipass==2.and.ISLAVE)  lni         = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = ginde_par(4,kfl_desti_par)
           if(ipass==2.and.ISLAVE)  lnb         = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = lneig_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  nneig       = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = slfbo_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  slfbo       = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = npoin_total
           if(ipass==2.and.ISLAVE)  npoin_total = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = nelem_total
           if(ipass==2.and.ISLAVE)  nelem_total = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = npoin_total
           if(ipass==2.and.ISLAVE)  npoin_total = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = nboun_total
           if(ipass==2.and.ISLAVE)  nboun_total = parin(ii)
           if(ipass==1) then
              npari = ii
              nparr = ir
              nparc = ic
              call memory_alloca(mem_servi(1:2,servi),'PARIN','par_sendat',parin,npari)
              call memory_alloca(mem_servi(1:2,servi),'PARRE','par_sendat',parre,nparr)
              if(ISLAVE) call par_receiv()
           end if
        end do
        if(IMASTER) call par_sendin()

        call memory_deallo(mem_servi(1:2,servi),'PARIN','par_sendat',parin)
        call memory_deallo(mem_servi(1:2,servi),'PARRE','par_sendat',parre)

     end if

  case(6_ip)
     !
     ! Exchange postprocess structure
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1

     do modul = 1,mmodu
        if( kfl_modul(modul) == 1 ) then
           do parii = 1,2
              npari = 0
              nparr = 0
              nparc = 0
              !
              ! Exchange of nsi_reaphy variables
              !
              call posdef(1_ip,dummi)
              !
              ! Allocate memory for the first pass
              !
              if(parii==1) then
                 call memory_alloca(mem_servi(1:2,servi),'PARIN','par_sendat',parin,npari)
                 call memory_alloca(mem_servi(1:2,servi),'PARRE','par_sendat',parre,nparr)
                 if( ISLAVE .or. kfl_ptask == 2 ) call par_broadc()
              end if
           end do

           if( IMASTER .and. kfl_ptask /= 2 ) call par_broadc()

           call memory_deallo(mem_servi(1:2,servi),'PARIN','par_sendat',parin)
           call memory_deallo(mem_servi(1:2,servi),'PARRE','par_sendat',parre)

        end if
     end do

     kfl_ptask = kfl_ptask_old

  case(7_ip)
     !
     ! Boundary conditions
     !
     if( IPARALL ) then

        strre='cderda_reaset'
        strin='cderda_reaset'
        strch='cderda_reaset'
        do parii=1,2
           npari=0
           nparr=0
           nparc=0
           !
           ! Other variables
           !
           call iexcha(kfl_schur)
           call iexcha(kfl_aiipr)

           if(parii==1) then
              call memory_alloca(mem_servi(1:2,servi),'PARIN','par_sendat',parin,npari)
              call memory_alloca(mem_servi(1:2,servi),'PARRE','par_sendat',parre,nparr)
              if( ISLAVE .or. kfl_ptask == 2 ) call par_broadc()
           end if
        end do

        if( IMASTER .and. kfl_ptask /= 2 ) call par_broadc()

        call memory_deallo(mem_servi(1:2,servi),'PARIN','par_sendat',parin)
        call memory_deallo(mem_servi(1:2,servi),'PARRE','par_sendat',parre)

     end if

  case(8_ip)

    if( IPARALL ) then

        strre='readim_reastr_reageo'
        strin='readim_reastr_reageo'
        strch='readim_reastr_reageo'
        do parii=1,2
           ipass = parii
           npari=0
           nparr=0
           nparc=0
           !
           ! Read in readim
           !
           npari=npari+1
           if(ipass==2.and.IMASTER) parin(npari) = npoin_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  npoin        = parin(npari)
           npari=npari+1
           if(ipass==2.and.IMASTER) parin(npari) = nelem_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  nelem        = parin(npari)
           npari=npari+1
           if(ipass==2.and.IMASTER) parin(npari) = nboun_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  nboun        = parin(npari)

           if( ipass == 1 ) then
              call memory_alloca(mem_servi(1:2,servi),'PARIN','par_sendat',parin,npari)
              call memory_alloca(mem_servi(1:2,servi),'PARRE','par_sendat',parre,nparr)
              if( ISLAVE ) call par_receiv()
           end if
        end do

        if( IMASTER ) call par_sendin()

        call memory_deallo(mem_servi(1:2,servi),'PARIN','par_sendat',parin)
        call memory_deallo(mem_servi(1:2,servi),'PARRE','par_sendat',parre)

     end if

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine par_sendat

