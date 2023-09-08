!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_sendat.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Exchange data
!> @details Exchange data
!> @} 
!------------------------------------------------------------------------
subroutine por_sendat(order)
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_porous
  use def_inpout
  use mod_memchk
  use mod_opebcs
  use mod_memory
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: kfl_ptask_old,dummi,ir,jr,kr,iwell

  nullify(parin)    ! This is needed here because someone has been untidy (not nullified) somewhere else
  nullify(parre)

  select case (order)

  case(1_ip)    

     !------------------------------------------------------------------- 
     !
     ! Exchange data read in por_reaphy, por_reanut and por_reaous
     !
     !------------------------------------------------------------------- 

     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     call Parall(29_ip)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of por_reaphy variables 
        !
        do kr=1,nprsa_por
           call iexcha(kfl_timei_por(kr))
        end do
        do kr=1,nprsa_por
           call iexcha(kfl_advec_por(kr))
        end do
        call iexcha(nmate_por)
        call iexcha(nwell_por)
        call iexcha(mrows_por)
        call iexcha(mroww_por)
        call iexcha(mheiw_por)

        call rexcha(comro_por)   
        call rexcha(comwa_por)   
        call rexcha(comoi_por)  
        call rexcha(bwref_por)   
        call rexcha(boref_por)   
        call rexcha(prref_por)   
        call rexcha(muwat_por)   
        call rexcha(muoil_por)   
        call rexcha(denwa_por)   
        call rexcha(denoi_por)   
        call rexcha(denhy_por)   
        do ir = 1,3
           call rexcha(gravi_por(ir))
        end do
        call rexcha(grnor_por)   
        call rexcha(prini_por)   
        !
        ! Exchange of por_reanut variables 
        !        
        call iexcha(kfl_ellen_por)
        do kr=1,nprsa_por
           call iexcha(kfl_sgsti_por(kr))
        end do
        do kr=1,nprsa_por
           call iexcha(kfl_sgsno_por(kr))
        end do
        do kr=1,nprsa_por
           call iexcha(kfl_taust_por(kr))
        end do
        do kr=1,nprsa_por
           call iexcha(kfl_ortho_por(kr))
        end do
        do kr=1,nprsa_por
           call iexcha(kfl_limit_por(kr))
        end do
        do kr=1,nprsa_por
           call iexcha(kfl_shock_por(kr))
        end do
        call iexcha(kfl_difwe_por)
        do kr=1,nprsa_por
           call iexcha(kfl_tiacc_por(kr))
        end do
        call iexcha(neule_por)
        call iexcha(kfl_tisch_por)
        call iexcha(kfl_normc_por)
        call iexcha(miinn_por)

        call rexcha(staco_por(1))
        call rexcha(staco_por(2))
        call rexcha(staco_por(3))
        call rexcha(shock_por)
        call rexcha(difwe_por)
        call rexcha(safet_por)
        call rexcha(sstol_por)
        call rexcha(cotol_por)
        call rexcha(relax_por)
        call rexcha(bemol_por)
        call rexcha(relsg_por) 
        call rexcha(tosgs_por) 
        solve_sol => solve(1:nprsa_por)
        call soldef(1_ip)
        !
        ! Exchange data read in por_reaous
        !
        call posdef(1_ip,dummi)
        !
        ! Allocate memory for the first pass
        !
        if( parii == 1 ) then
           call memory_alloca(mem_servi(1:2,ID_PARALL),'PARIN','por_sendat',parin,npari)
           call memory_alloca(mem_servi(1:2,ID_PARALL),'PARRE','por_sendat',parre,nparr)
           if( ISLAVE .or. kfl_ptask == 2 ) call Parall(2_ip)
        end if
     end do
     if( IMASTER .and. kfl_ptask /= 2 ) call Parall(2_ip)
     call memory_deallo(mem_servi(1:2,ID_PARALL),'PARIN','por_sendat',parin )
     call memory_deallo(mem_servi(1:2,ID_PARALL),'PARRE','por_sendat',parre )
     !
     ! Allocatable arrays
     !
     call Parall(30_ip)
     if( ISLAVE ) then
        call por_memphy(1_ip)   ! For the moment does nothing
        call por_memphy(3_ip)   ! TKREL_POR & NROWS_POR
        call por_memphy(4_ip)   ! well
     end if
     !
     ! Physical properties
     !
     do parii = 1,2 
        npari = 0
        nparr = 0
        nparc = 0
        do iwell = 1,nwell_por
           call iexcha(tywel_por(iwell) % itype)
           call iexcha(tywel_por(iwell) % ntime)
           call iexcha(tywel_por(iwell) % nposi)
           call rexcha(tywel_por(iwell) % radiu)
        end do
        if( parii == 1 ) then
           call memory_alloca(mem_servi(1:2,ID_PARALL),'PARIN','por_sendat',parin,npari)
           call memory_alloca(mem_servi(1:2,ID_PARALL),'PARRE','por_sendat',parre,nparr)
           if( ISLAVE .or. kfl_ptask == 2 ) call Parall(2_ip)
        end if
     end do
     if( IMASTER .and. kfl_ptask /= 2 ) call Parall(2_ip)
     call memory_deallo(mem_servi(1:2,ID_PARALL),'PARIN','por_sendat',parin )
     call memory_deallo(mem_servi(1:2,ID_PARALL),'PARRE','por_sendat',parre )     
     !
     ! Allocate wells
     !
     if( ISLAVE ) then
        do iwell = 1,nwell_por
           igene = iwell
           call por_memphy(5_ip)
        end do
     end if
     !
     ! Physical properties
     !
     do parii = 1,2 
        npari = 0
        nparr = 0
        nparc = 0
        !
        ! Exchange of por_reaphy variables whose dimensions depend
        ! on what is read in por_reaphy
        !
        do iwell = 1,nwell_por
           do ir = 1,tywel_por(iwell) % ntime
              do jr = 1,2
                 call rexcha(tywel_por(iwell) % q_table(ir,jr))
              end do
              do jr = 1,tywel_por(iwell) % nposi + 1
                 call rexcha(tywel_por(iwell) % pbh_table(ir,jr))
              end do
           end do
           do ir = 1,tywel_por(iwell) % nposi
              call rexcha(tywel_por(iwell) % pbh_coord(ir))
           end do
        end do

        do kr=1,nmate_por
           call iexcha(nrows_por(kr))
        end do
        do kr=1,nwell_por
           call iexcha(nroww_por(kr))
        end do
        do kr=1,nwell_por
           call iexcha(nheiw_por(kr))
        end do
        do kr=1,nwell_por
           call iexcha(kfl_wellc_por(kr))
        end do
        do ir=1,nmate_por
           do jr=1,mrows_por
              do kr=1,3
                 call rexcha(tkrel_por(kr,jr,ir))
              end do
           end do
        end do
        do ir=1,nwell_por
           do jr=1,mroww_por
              do kr=1,mheiw_por+1
                 call rexcha(wvalu_por(kr,jr,ir))
              end do
           end do
        end do
        do kr=1,nwell_por
           call rexcha(rwell_por(kr))
        end do

        if( parii == 1 ) then
           call memory_alloca(mem_servi(1:2,ID_PARALL),'PARIN','por_sendat',parin,npari)
           call memory_alloca(mem_servi(1:2,ID_PARALL),'PARRE','por_sendat',parre,nparr)
           if( ISLAVE .or. kfl_ptask == 2 ) call Parall(2_ip)
        end if
     end do
     if( IMASTER .and. kfl_ptask /= 2 ) call Parall(2_ip)
     call memory_deallo(mem_servi(1:2,ID_PARALL),'PARIN','por_sendat',parin )
     call memory_deallo(mem_servi(1:2,ID_PARALL),'PARRE','por_sendat',parre )     
     !
     ! Deallocate wells for Master
     !
     if( IMASTER ) then
        call por_memphy(-5_ip)
     end if

     !------------------------------------------------------------------- 
     !
     ! Variables read in reabcs
     !
     !------------------------------------------------------------------- 

     call Parall(27_ip)
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of por_reabcs variables 
        !
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           call memory_alloca(mem_servi(1:2,ID_PARALL),'PARIN','por_sendat',parin,npari)
           call memory_alloca(mem_servi(1:2,ID_PARALL),'PARRE','por_sendat',parre,nparr)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(2_ip)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(2_ip)
     call memory_deallo(mem_servi(1:2,ID_PARALL),'PARIN','por_sendat',parin )
     call memory_deallo(mem_servi(1:2,ID_PARALL),'PARRE','por_sendat',parre )       
     !
     ! Boundary codes
     !
     call spnbcs(tncod_por)

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine por_sendat
