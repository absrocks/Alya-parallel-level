subroutine wav_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_sendat
  ! NAME
  !    wav_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    wav_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_wavequ
  use def_inpout
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,jc,ki,kr,kc
  integer(ip)             :: ielem,ipoin,dummi
  integer(ip)             :: jnode,jpoin,jelem,inode
  integer(ip)             :: ibcas,ixchn,kfl_ptask_old 
  integer(4)              :: istat

  ibcas= 2_ip
  ixchn= 300_ip


  select case (order)

  case(1)     
     !
     ! Exchange data read in wav_reaphy, wav_reanut and wav_reaous
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     call Parall(27_ip)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of wav_reaphy variables 
        !
        call iexcha(kfl_sourc_wav)
        do jr=1,nsour_wav
           call rexcha(sourc_wav(jr)) 
        end do        
        call iexcha(nmate_wav)
        !
        ! Exchange of wav_reanut variables 
        !        
        call iexcha(kfl_tiacc_wav)
        call iexcha(kfl_tisch_wav)
        call iexcha(kfl_timet_wav)
        call iexcha(kfl_subgs_wav)
        call iexcha(kfl_massm_wav)
        call iexcha(neule_wav)
        call rexcha(safet_wav) 
        call rexcha(sstol_wav) 
        call rexcha(nebet_wav)
        call rexcha(negam_wav)
        solve_sol => solve
        call soldef(1_ip)
        !
        ! Exchange data read in wav_reaous
        !
        call posdef(1_ip,dummi)
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','wav_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','wav_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(ibcas)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(ibcas)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','wav_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','wav_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','wav_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','wav_sendat',0_ip)     
     !
     ! Allocatable arrays: LMATE_WAV and LMATN_WAV
     !
     if(kfl_paral>=1) call wav_memphy()
     party =  1 ! vector dimensioned nelem
     pardi =  1 ! = 1, 2, 3 -->  number of columns of a vector (pari1, 2,3 )
     parki =  1 ! = 1, 2, 3 -->  integer, real or character
     pari1 => lmate_wav        
     call Parall(ixchn)

     kfl_ptask = kfl_ptask_old
     !
     ! Physical properties
     !
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of wav_reaphy variables 
        !
        do ji=1,nmate_wav
           call iexcha(lawde_wav(ji)) 
        end do    
        do ji=1,nmate_wav
           call iexcha(lawka_wav(ji)) 
        end do    
        do jr=1,nmate_wav
           do kr=1,ncoef_wav
              call rexcha(densi_wav(kr,jr))
           end do
        end do    
        do jr=1,nmate_wav
           do kr=1,ncoef_wav
              call rexcha(kappa_wav(kr,jr))
           end do
        end do            
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','wav_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','wav_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(ibcas)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(ibcas)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','wav_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','wav_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','wav_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','wav_sendat',0_ip)     
 
  case(2)

     call Parall(29_ip)
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of wav_reabcs variables 
        !
        call iexcha(kfl_onnod_wav)
        call iexcha(kfl_onbou_wav)    
        call iexcha(kfl_absor_wav)    
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','wav_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','wav_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(ibcas)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(ibcas)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','wav_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','wav_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','wav_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','wav_sendat',0_ip)     

     if(kfl_onnod_wav==1) then
        !
        ! Boundary conditions on nodes
        !
        if(kfl_paral>=1) call wav_membcs(1_ip)
        
        if(kfl_paral/=0.or.kfl_ptask/=2) then
           !
           ! KFL_FIXNO_WAV
           !
           party =  3 ! node
           pardi =  1 ! # dimensions
           parki =  1 ! kind=integer
           pari1 => kfl_fixno_wav
           call Parall(ixchn)
           !
           ! KFL_BVESS_WAV
           !
           party =  3 ! node
           pardi =  1 ! # dimensions
           parki =  2 ! kind=real
           parr1 => bvess_wav
           call Parall(ixchn)
           
        end if

     end if

     if(kfl_onbou_wav==1) then
        !
        ! Boundary conditions on boundaries
        !
        if(kfl_paral>=1) call wav_membcs(2_ip)
        
        if(kfl_paral/=0.or.kfl_ptask/=2) then
           !
           ! KFL_FIXBO_WAV
           !
           party =  2 ! boundaries
           pardi =  1 ! # dimensions
           parki =  1 ! kind=integer
           pari1 => kfl_fixbo_wav
           call Parall(ixchn)
           
        end if

     end if

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine wav_sendat
