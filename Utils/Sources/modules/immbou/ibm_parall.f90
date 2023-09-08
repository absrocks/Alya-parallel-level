subroutine ibm_parall(order)
  !-----------------------------------------------------------------------
  !****f* Temper/ibm_parall
  ! NAME
  !    ibm_parall
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_inpout
  use mod_memchk
  use def_immbou
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,kfl_ptask_old,dummi
  integer(4)              :: istat

  if ( ISEQUEN ) return

  select case ( order )

  case ( 1_ip )    

     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1

     !-------------------------------------------------------------------
     !
     ! Exchange NIMBO and allocate memory for particles
     !
     !-------------------------------------------------------------------

     call parari('BCT',0_ip,1_ip,nimbo)
     if( ISLAVE ) call ibm_defini(0_ip)

     !-------------------------------------------------------------------
     !
     ! Dimensions and and solver
     !
     !-------------------------------------------------------------------

     strre = 'ibm_parall'
     strin = 'ibm_parall'
     strch = 'ibm_parall'
     do parii = 1,2 
        npari = 0
        nparr = 0
        nparc = 0
        !
        ! Variables read ibm_reageo, ibm_reaphy and ibm_reanut
        !
        call ibm_defini(1_ip)
        !
        ! Variables computed in ibm_readim and ibm_cderda
        !
        call iexcha(ndime)
        do ji=1,nelty
           call iexcha(lexib(ji))
        end do
        do ji=1,nelty
           call iexcha(ngaib(ji))
        end do
        do ji=1,nelty
           call iexcha(lruib(ji))
        end do
        do ji=1,nelty
           call iexcha(lquib(ji))
        end do
        call iexcha(mnoib)
        call iexcha(mnodi)
        call iexcha(mgaib)
        call iexcha(mgaui)
        !
        ! Solver
        !
        solve_sol => solve(1:1)
        call soldef(1_ip)
        !
        ! Postprocess
        !
        call posdef(1_ip,dummi)

        if( parii == 1 ) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'PARIN','sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'PARRE','sendat',parre)
           if( ISLAVE ) call Parall(two)
        end if

     end do

     if( IMASTER ) call Parall(two)

     call memchk(two,istat,mem_modul(1:2,modul),'PARIN','sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'PARIN','sendat',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'PARRE','sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'PARRE','sendat',0_ip)

     !-------------------------------------------------------------------
     !
     ! COOIN, COOIB, COOI2, LNOIB, LTYIB 
     !
     !-------------------------------------------------------------------

     if( ISLAVE ) call ibm_memall(1_ip)

     do parii = 1,2 
        npari = 0
        nparr = 0
        nparc = 0

        call ibm_defini(7_ip)

        if( parii == 1 ) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'PARIN','ibm_parall',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'PARRE','ibm_parall',parre)
           if( ISLAVE ) call Parall(two)
        end if

     end do

     if( IMASTER ) call Parall(two)

     npari = 0
     nparr = 0
     nparc = 0
     call memchk(two,istat,mem_modul(1:2,modul),'PARIN','ibm_parall',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'PARIN','ibm_parall',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'PARRE','ibm_parall',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'PARRE','ibm_parall',0_ip)

     kfl_ptask = kfl_ptask_old

  end select

end subroutine ibm_parall

