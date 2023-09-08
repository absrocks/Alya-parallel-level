subroutine pts_parall(order)
  !-----------------------------------------------------------------------
  !****f* Parall/pts_parall
  ! NAME
  !    pts_parall
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
  use def_partis
  use mod_opebcs
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_BROADCAST
  use mod_communications, only : PAR_EXCHANGE
  use mod_pts_injection,  only : pts_injection_parallelization
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ii,dummi,jj

  if( ISEQUEN ) return

  select case ( order )

  case( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Broadcast data read in *.pts.dat file
     !
     !-------------------------------------------------------------------

     do parii = 1,2 
        npari = 0 
        nparr = 0
        nparc = 0
        nparl = 0
        !
        ! Physical problem
        !
        call PAR_EXCHANGE(mlagr,             parin,npari,parii)
        call PAR_EXCHANGE(dtmin_pts,         parre,nparr,parii)
        call PAR_EXCHANGE(dimin_pts,         parre,nparr,parii)
        call PAR_EXCHANGE(mean_free_path_pts,parre,nparr,parii)           
        !
        ! Type description
        !
        do ii = 1,mtyla
           call PAR_EXCHANGE(parttyp(ii) % kfl_exist,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_modla,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_grafo,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_buofo,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_drafo,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_extfo,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_brown,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_saffm,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_schem,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % denpa,    parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % spher,    parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % diame,    parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % calor,    parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % emisi,    parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % scatt,    parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % diffu,    parre,nparr,parii) 
           call PAR_EXCHANGE(mlapr,parttyp(ii) % prope,parre,nparr,parii) 
           call PAR_EXCHANGE(mlapr,parttyp(ii) % prova,parin,npari,parii) 
        end do
        !
        ! Numerical problem
        !
        call PAR_EXCHANGE(kfl_adapt_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_usbin_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_order_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_walld_pts,parin,npari,parii)
        call PAR_EXCHANGE(gamma_pts,parre,nparr,parii)
        call PAR_EXCHANGE(beta_pts ,parre,nparr,parii)
        call PAR_EXCHANGE(chale_pts,parre,nparr,parii)
        call PAR_EXCHANGE(safet_pts,parre,nparr,parii)
        call PAR_EXCHANGE(dtime_pts,parre,nparr,parii)
        solve_sol => solve(1:)
        call soldef(1_ip)
        !
        ! Output and postprocess
        !                  
        call posdef(1_ip,dummi)

        call PAR_EXCHANGE(kfl_ptsres_binary        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_ptsres_split         ,parin,npari,parii)


        call PAR_EXCHANGE(kfl_posla_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_oudep_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_depos_surface_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_oufre_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_dbfre_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_exacs_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_rstar_pts        ,parin,npari,parii)

#ifdef DBPARTICLES           
        call PAR_EXCHANGE(dbSettings % kfl_db_deep      ,parin,npari,parii)
        call PAR_EXCHANGE(dbSettings % maxInserts       ,parin,npari,parii)
        call PAR_EXCHANGE(dbSettings % port             ,parin,npari,parii)
        call PAR_EXCHANGE(dbSettings % numParticlesBatch,parin,npari,parii)
        call PAR_EXCHANGE(dbSettings % kfl_db_url_conn  ,parch,nparc,parii)           
#endif
        call PAR_EXCHANGE(mvarp_pts,postprocess_var_pts,parlo,nparl,parii)
        call PAR_EXCHANGE(mvarp_pts,deposition_var_pts ,parlo,nparl,parii)
        !
        ! Boundary conditions
        !
        call PAR_EXCHANGE(veloc_field_id, parin,npari,parii)
        call PAR_EXCHANGE(pts_minj,kfl_injla_pts,parin,npari,parii)
        call PAR_EXCHANGE(pts_minj,injector_particle_distribution,parin,npari,parii)
        call PAR_EXCHANGE(pts_minj,injector_npts_asis,parin,npari,parii)
        call PAR_EXCHANGE(pts_minj,kfl_injty_pts ,parin,npari,parii)
        call PAR_EXCHANGE(pts_minj,kfl_random_pts,parin,npari,parii)              
        call PAR_EXCHANGE(kfl_boundary_injection,parin,npari,parii) 
        call PAR_EXCHANGE(pts_minj,codbo_pts,parin,npari,parii) 
        call PAR_EXCHANGE(kfl_injve_pts,parin,npari,parii) 
        call PAR_EXCHANGE(tinla_pts,parre,nparr,parii)
        call PAR_EXCHANGE(tfila_pts,parre,nparr,parii)
        call PAR_EXCHANGE(tpela_pts,parre,nparr,parii)
        call PAR_EXCHANGE(mpala,parla2_pts,parre,nparr,parii)

        do ii = 1,pts_minj
           do jj = 1,mpala
              call PAR_EXCHANGE(parla_pts(ii,jj),parre,nparr,parii)
           end do
        end do

        if( parii == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'PARIN','pts_parall',parin,npari)
           call memory_alloca(mem_modul(1:2,modul),'PARRE','pts_parall',parre,nparr)
           call memory_alloca(mem_modul(1:2,modul),'PARLO','pts_parall',parlo,nparl)
           if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
           if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
           if( ISLAVE  ) call PAR_BROADCAST(parlo,      'IN MY CODE')
           if( ISLAVE  ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
        else
           if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
           if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
           if( IMASTER ) call PAR_BROADCAST(parlo,      'IN MY CODE')
           if( IMASTER ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
        end if

     end do

     call memory_deallo(mem_modul(1:2,modul),'PARIN','pts_parall',parin)
     call memory_deallo(mem_modul(1:2,modul),'PARRE','pts_parall',parre)
     call memory_deallo(mem_modul(1:2,modul),'PARLO','pts_parall',parlo)

  end select

  npari = 0
  nparr = 0
  nparc = 0
  nparl = 0
  !
  ! Injector with coordinates
  !
  call pts_injection_parallelization()
  !
  ! Boundary conditions
  !
  call spbbcs(tbcod_pts)

end subroutine pts_parall

