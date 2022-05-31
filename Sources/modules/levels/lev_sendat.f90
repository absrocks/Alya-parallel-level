subroutine lev_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_sendat
  ! NAME
  !    lev_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    lev_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master 
  use def_kermod
  use def_solver
  use def_domain
  use def_levels
  use def_inpout
  use mod_memchk
  use mod_opebcs
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,ki,dummi
  integer(ip)             :: ixchn,kfl_ptask_old 
  integer(4)              :: istat

  select case (order)

  case(1_ip)     
     !
     ! Exchange data read in lev_reaphy, lev_reanut and lev_reaous
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of lev_reaphy variables 
        !
        call iexcha(kfl_inlev_lev)
        call iexcha(kfl_advec_lev)
        call iexcha(kfl_reave_lev)
        call iexcha(nmate_lev)
        call rexcha(thicl)                    ! belongs to defmaster. Shared with other modules, but only after itask_turnon.
        !
        ! Exchange of lev_reanut variables 
        !        
        call iexcha(kfl_timet_lev)
        call iexcha(kfl_tisch_lev)
        call iexcha(kfl_timco_lev)
        call iexcha(kfl_tiacc_lev)
        call iexcha(kfl_normc_lev)
        call iexcha(kfl_ellen_lev)
        call iexcha(kfl_zonal_lev)
        call iexcha(neule_lev)
        call iexcha(miinn_lev)
        call iexcha(inred_lev)
        call iexcha(nfred_lev)
        call iexcha(tyred_lev)
        call iexcha(nstre_lev)
        call iexcha(kfl_locre_lev)
        call iexcha(nbitr_lev)
        call iexcha(kfl_corvo_lev)
        call rexcha(safet_lev) 
        call rexcha(sstol_lev)
        call rexcha(cotol_lev)
        call rexcha(cpuit_lev)
        call rexcha(supgp_lev)
        solve_sol => solve
        call soldef(1_ip)
        !
        ! Exchange data read in lev_reaous
        !
        call posdef(1_ip,dummi)
        call iexcha(npp_gauge_lev)           
        call iexcha(npp_nbgau_lev)
        call iexcha(npp_inter_lev)
        do ji=1,ngaug_lev
           call iexcha(typga_lev(ji)) 
        end do
        do ki=1,ngaug_lev
           do ji=1,3
              call rexcha(cogau_lev(ji,ki))
           end do
        end do
        !
        ! Exchange of lev_reabcs variables 
        !
        call iexcha(kfl_inlev_lev) 
        call iexcha(kfl_conbc_lev) 
        call rexcha(height_lev)
        do ji=1, 4
           call rexcha(lev_dam_break(ji))
        end do
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','lev_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','lev_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call par_broadc()
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'parin','lev_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','lev_sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'parre','lev_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','lev_sendat',0_ip)     
 
     !------------------------------------------------------------------- 
     !
     ! Variables read in reabcs
     !
     !------------------------------------------------------------------- 

     call spnbcs(tncod_lev)
     call spgbcs(tgcod_lev)

  case (6)
        
      !  call par_reduce()  
!!$        call lev_memall(2_ip)                   
      !  call par_redvec()                             

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine lev_sendat
