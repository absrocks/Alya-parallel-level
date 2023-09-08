subroutine rad_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_sendat
  ! NAME
  !    rad_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    rad_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_radiat
  use def_inpout
  use mod_memchk
  use mod_opebcs
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,ki,kr,ir,kfl_ptask_old,dummi
  integer(4)              :: istat

  select case (order)

  case(1_ip)    

     !------------------------------------------------------------------- 
     !
     ! Exchange data read in rad_reaphy, rad_reanut and rad_reaous
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
        ! Exchange of rad_reaphy variables 
        !
         call iexcha(kfl_atest_rad)
         call iexcha(kfl_parti_rad)
         call iexcha(idtem_rad)

        !
        ! Exchange of rad_reanut variables 
        !        

        call iexcha(kfl_ellen_rad)
        call iexcha(kfl_sgsno_rad)
        call iexcha(kfl_taust_rad)
        call iexcha(kfl_ortho_rad)
        call iexcha(kfl_limit_rad)

        call iexcha(kfl_bubbl_rad)
        call iexcha(kfl_assem_rad)
        call iexcha(kfl_normc_rad)
        call iexcha(miinn_rad)

        call rexcha(staco_rad(1))
        call rexcha(staco_rad(2))
        call rexcha(staco_rad(3))
        call rexcha(cotol_rad)
        call rexcha(relsg_rad)
        call rexcha(tosgs_rad) 
        solve_sol => solve(1:1)
        call soldef(1_ip)
        !
        ! Exchange data read in rad_reaous
        !
        call posdef(1_ip,dummi)
        call iexcha(kfl_splot_rad)
        call iexcha(kfl_psmat_rad)
        call iexcha(kfl_exacs_rad)
        call iexcha(npp_bound_rad)
        do ji=1,nexap_rad
           call rexcha(expar_rad(ji))
        end do
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','rad_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','rad_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(2_ip)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','rad_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','rad_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','rad_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','rad_sendat',0_ip)     
     !
     ! Allocatable arrays: Scatt_rad and absor_rad and aniso_rad
     !
     call Parall(30_ip)
     if(kfl_paral>=1) then
        call rad_memphy(1_ip)
        call rad_memphy(2_ip)
     end if

     kfl_ptask = kfl_ptask_old
     !
     ! Physical properties
     !
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of rad_reaphy variables whose dimensions depend
        ! on what is read in rad_reaphy
        !
        do ji=1,nspec_rad
           call rexcha(scatt_rad(ji)) 
        end do
        do ji=1,nspec_rad
           call rexcha(absor_rad(ji)) 
        end do
        do ji=1,nspec_rad
           call rexcha(aniso_rad(ji)) 
        end do

        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','rad_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','rad_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(2_ip)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','rad_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','rad_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','rad_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','rad_sendat',0_ip)     

     call Parall(27_ip)
 
     !------------------------------------------------------------------- 
     !
     ! Variables read in reabcs
     !
     !------------------------------------------------------------------- 

     call spnbcs(tncod_rad)
     !call spnbcs(tgcod_rad)
     call spbbcs(tbcod_rad)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of rad_reabcs variables 
        !
        call iexcha(kfl_conbc_rad) 
        call iexcha(kfl_inico_rad)
        call iexcha(kfl_intbc_rad)
        call iexcha(npnat_rad)
        do jr=1,nmate
           call rexcha(bvcoe_rad(jr))
        end do
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','rad_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','rad_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(2_ip)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','rad_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','rad_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','rad_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','rad_sendat',0_ip)     

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine rad_sendat
