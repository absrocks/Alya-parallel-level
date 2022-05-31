subroutine ale_parall(order)
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_parall
  ! NAME
  !    ale_parall
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    ale_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_alefor
  use def_inpout
  use mod_memchk
  use mod_opebcs
  use def_kermod,         only : number_space_time_function
  use mod_alefor,         only : alefor_memory_allocate
  use mod_communications, only : PAR_BROADCAST
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: kfl_ptask_old,dummi,ifunc,ji
  integer(4)              :: istat

  if( ISEQUEN ) return
  
  select case (order)

  case(1_ip)    

     !------------------------------------------------------------------- 
     !
     ! Exchange data read in ale_reaphy, ale_reanut and ale_reaous
     !
     !------------------------------------------------------------------- 

     call PAR_BROADCAST(nrbod)
     call PAR_BROADCAST(kfl_rigid_ale)
     if( ISLAVE .and. kfl_rigid_ale == 1 ) call alefor_memory_allocate('RBBOU')
        
     strre = 'ale_parall'
     strin = 'ale_parall'
     strch = 'ale_parall'
     
     do parii=1,2 
        npari=0;nparr=0;nparc=0
        !
        ! Exchange of ale_reaphy, ale_reanut variables 
        !
        call iexcha(kfl_smoot_ale)
        call iexcha(kfl_timef_ale)
        call iexcha(nsmoo_ale)
        call iexcha(kfl_defor_ale)
        call iexcha(ndefo_ale)
        call iexcha(kfl_smobo_ale)
        call iexcha(kfl_fixsm_ale)
        call iexcha(nsmob_ale)
        call iexcha(kfl_crist_ale)
        call iexcha(kfl_foexo_ale)
        call iexcha(kfl_disor_ale)
        call iexcha(kfl_nforc_ale)
        call iexcha(kfl_sensi_ale)
        call iexcha(moddi_ale)
        call iexcha(modvi_ale)
        call rexcha(ansmo_ale)
        call rexcha(resmo_ale)
        call rexcha(defor_param_ale)
        !
        ! Variables for rigid body
        !
        if( kfl_rigid_ale == 1 ) call ale_sendat(1_ip)
        !
        ! Variables computed in ibm_readim and ibm_cderda
        !
        do ji=1,nelty
           call iexcha(lexib(ji))
        end do
        do ji=1,nelty
           call iexcha(ngaib(Ji))
        end do
        do ji=1,nelty
           call iexcha(lruib(ji))
        end do

        call iexcha(mnoib)
        call iexcha(mnodi)
        call iexcha(mgaib)
        !
        ! Solver
        !
        solve_sol => solve(1:1)
        call soldef(1_ip)
        !
        ! Exchange data read in ale_reaous
        !
        call posdef(1_ip,dummi)
        !
        ! Allocate memory for the first pass
        !
        if( parii == 1 ) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parin','ale_parall',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'parre','ale_parall',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call par_broadc()
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call par_broadc()

     call memchk(two,istat,mem_modul(1:2,modul),'PARIN','ale_parall',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'PARIN','ale_parall',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'PARRE','ale_parall',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'PARRE','ale_parall',0_ip)

     call spnbcs(tncod_ale)
     call spgbcs(tgcod_ale)
     call spbbcs(tbcod_ale)

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine ale_parall
