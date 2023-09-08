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
  use def_kermod, only: number_space_time_function
  use mod_memchk
  use mod_opebcs
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: kfl_ptask_old,dummi,ifunc,ji
  integer(4)              :: istat

  select case (order)

  case(1_ip)    

     !------------------------------------------------------------------- 
     !
     ! Exchange data read in ale_reaphy, ale_reanut and ale_reaous
     !
     !------------------------------------------------------------------- 

     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     call Parall(29_ip)

     !-------------------------------------------------------------------
     !
     ! Exchange nrbod and allocate memory for particles
     !
     !-------------------------------------------------------------------

     call parari('BCT',0_ip,1_ip,nrbod)
     call parari('BCT',0_ip,1_ip,kfl_rigid_ale)
     if( ISLAVE .and. kfl_rigid_ale == 1 ) call ale_sendat(0_ip)

     strre = 'ale_parall'
     strin = 'ale_parall'
     strch = 'ale_parall'
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
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
        call iexcha(kfl_suppo_ale)
        call iexcha( kfl_sensi_ale )
        call iexcha( moddi_ale )
        call iexcha( modvi_ale )
        call rexcha(ansmo_ale)
        call rexcha(resmo_ale)
        !
        ! Variables for rigid body
        !
        if( kfl_rigid_ale == 1 ) call ale_sendat(1_ip)
        !
        !
        ! Variables computed in ibm_readim and ibm_cderda
        !
#ifndef NDIMEPAR
        call iexcha(ndime)
#endif
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
        ! Exchange of ale_reabcs variables 
        !
        call iexcha(kfl_conbc_ale)
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
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','ale_parall',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','ale_parall',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(2_ip)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'PARIN','ale_parall',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'PARIN','ale_parall',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'PARRE','ale_parall',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'PARRE','ale_parall',0_ip)

     if ( kfl_rigid_ale == 1 ) then 
        !-------------------------------------------------------------------
        !
        ! NPOIB, NBOIB, MASSA etc 
        !
        !-------------------------------------------------------------------

        do parii = 1,2 
           npari = 0
           nparr = 0
           nparc = 0

           call ale_sendat(2_ip)
           call ale_sendat(3_ip)    ! only those read in restart

           if( parii == 1 ) then
              allocate(parin(npari),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'PARIN','ale_parall',parin)
              allocate(parre(nparr),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'PARRE','ale_parall',parre)
              if( ISLAVE ) call Parall(two)
           end if

        end do

        if( IMASTER ) call Parall(two)

        npari = 0
        nparr = 0
        nparc = 0
        call memchk(two,istat,mem_modul(1:2,modul),'PARIN','ale_parall',parin)
        deallocate(parin,stat=istat)
        if(istat/=0) call memerr(two,'PARIN','ale_parall',0_ip)
        call memchk(two,istat,mem_modul(1:2,modul),'PARRE','ale_parall',parre)
        deallocate(parre,stat=istat)
        if(istat/=0) call memerr(two,'PARRE','ale_parall',0_ip)

        !-------------------------------------------------------------------
        !
        ! COOIN, COOIB, COOI2, LNOIB, LTYIB 
        !
        !-------------------------------------------------------------------

        do parii = 1,2 
           npari = 0
           nparr = 0
           nparc = 0

           call ale_sendat(7_ip)

           if( parii == 1 ) then
              allocate(parin(npari),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'PARIN','ale_parall',parin)
              allocate(parre(nparr),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'PARRE','ale_parall',parre)
              if( ISLAVE ) call Parall(two)
           end if

        end do

        if( IMASTER ) call Parall(two)

        npari = 0
        nparr = 0
        nparc = 0
        call memchk(two,istat,mem_modul(1:2,modul),'PARIN','ale_parall',parin)
        deallocate(parin,stat=istat)
        if(istat/=0) call memerr(two,'PARIN','ale_parall',0_ip)
        call memchk(two,istat,mem_modul(1:2,modul),'PARRE','ale_parall',parre)
        deallocate(parre,stat=istat)
        if(istat/=0) call memerr(two,'PARRE','ale_parall',0_ip)

     end if

     call Parall(30_ip)
     kfl_ptask = kfl_ptask_old

     call spnbcs(tncod_ale)
     call spgbcs(tgcod_ale)
     call spbbcs(tbcod_ale)
     !
     ! KFL_FUNTY_ALE
     !
     if( kfl_conbc_ale == 0 .and. number_space_time_function==0_ip ) then
        dummi = 2 * mfunc_ale
        strin = 'kfl_funty_ale'
        if( ISLAVE )  call ale_membcs(47_ip)
        call parari('BCT',0_ip,dummi,kfl_funty_ale)
        !
        ! KFL_FUNPA_ALE
        !
        strre = 'funpa_ale'
        do ifunc = 1,mfunc_ale
           if( kfl_funty_ale(ifunc,1) /= 0 ) then
              igene = ifunc
              dummi = kfl_funty_ale(ifunc,2)
              if( ISLAVE )  call ale_membcs( 46_ip)
              call pararr('BCT',0_ip,dummi,funpa_ale(ifunc) % a)
              if( IMASTER ) call ale_membcs(-46_ip)
           end if
        end do
        if( IMASTER ) call ale_membcs(-47_ip)
     end if
     !
     ! Support geometry
     !
     if( kfl_suppo_ale == 1 ) then     
        call parari('BCT',0_ip,1_ip,npoin_ad)
        call parari('BCT',0_ip,1_ip,nboun_ad)
        call parari('BCT',0_ip,1_ip,mnodb_ad)
        if( ISLAVE )  call ale_membcs(5_ip)
        call pararr('BCT',0_ip,npoin_ad*ndime,   coord_ad)
        call parari('BCT',0_ip,nboun_ad*mnodb_ad,lnodb_ad)
        call parari('BCT',0_ip,nboun_ad         ,ltypb_ad)
        if( IMASTER ) call ale_membcs(-5_ip) 
     end if

  case(2_ip)    

     !------------------------------------------------------------------- 
     !
     ! Exchange data read in restart
     !
     !------------------------------------------------------------------- 

     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     call Parall(29_ip)     ! por ahora lo dejo idem al case 1

     strre = 'ale_parall'
     strin = 'ale_parall'
     strch = 'ale_parall'

     if ( kfl_rigid_ale == 1 ) then 

        !-------------------------------------------------------------------
        !
        ! Only those values read in restart
        !
        !-------------------------------------------------------------------

        do parii = 1,2 
           npari = 0
           nparr = 0
           nparc = 0

           call ale_sendat(3_ip)    ! only those read in restart

           if( parii == 1 ) then
              allocate(parin(npari),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'PARIN','ale_parall',parin)
              allocate(parre(nparr),stat=istat)
              call memchk(zero,istat,mem_modul(1:2,modul),'PARRE','ale_parall',parre)
              if( ISLAVE ) call Parall(two)
           end if

        end do

        if( IMASTER ) call Parall(two)

        npari = 0
        nparr = 0
        nparc = 0
        call memchk(two,istat,mem_modul(1:2,modul),'PARIN','ale_parall',parin)
        deallocate(parin,stat=istat)
        if(istat/=0) call memerr(two,'PARIN','ale_parall',0_ip)
        call memchk(two,istat,mem_modul(1:2,modul),'PARRE','ale_parall',parre)
        deallocate(parre,stat=istat)
        if(istat/=0) call memerr(two,'PARRE','ale_parall',0_ip)

     end if

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine ale_parall
