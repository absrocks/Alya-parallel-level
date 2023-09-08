subroutine qua_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_sendat
  ! NAME
  !    qua_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    qua_turnon
  !    qua_parallel
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_quanty
  use def_inpout
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,ki,kr,ir,kfl_ptask_old,ii
  integer(4)              :: istat

  select case (order)

  case(1_ip)     
     !
     ! Exchange data read in tem_reaphy, tem_reanut and tem_reaous
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     call Parall(29_ip)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of qua_reaphy variables 
        !
        call iexcha(kfl_timei_qua  )
        call iexcha(kfl_dftgs_qua  )
        call iexcha(kfl_potxc_qua  )
        call iexcha(kfl_alele_qua  )        
        call iexcha(kfl_nolin_qua  )
        call iexcha(kfl_perio_qua  )

        call iexcha(kfl_coulo_qua  )
        call iexcha(kfl_bfieldx_qua)
        call iexcha(kfl_bfieldy_qua)
        call iexcha(kfl_bfieldz_qua)
        call iexcha(klf_btemp_law  )
        call iexcha(kfl_efieldx_qua)
        call iexcha(kfl_efieldy_qua)
        call iexcha(kfl_efieldz_qua)
        call iexcha(klf_etemp_law  )

        call iexcha(kfl_spinb_qua  )
        call iexcha(kfl_spiorb_qua )
        call iexcha(kfl_relat_qua  )
        call iexcha(kfl_vother_qua )

        call iexcha(ncuanpal_qua   )
        call iexcha(lcuanorb_qua   )
        call iexcha(nspin_qua      )
        call iexcha(lawma_qua      )

        call iexcha(law_vother_qua)
        call iexcha(natoms_qua)
        call iexcha(nespecies_qua)

        call iexcha(ncomp_eig)
        call iexcha(nestates)
        call iexcha(noutput)    

        call rexcha(eig_evol_qua  )
        call rexcha(coulo_qua     )
        call rexcha(bfieldx_qua   )
        call rexcha(bfieldy_qua   )
        call rexcha(bfieldz_qua   )
        call rexcha(btemp_law_qua )
        call rexcha(efieldx_qua   )
        call rexcha(efieldy_qua   )
        call rexcha(efieldz_qua   )
        call rexcha(etemp_law_qua )
        call rexcha(frecuencie    )
        call rexcha(spinb_qua     )
        call rexcha(spiorb_qua    )
        call rexcha(relat_qua     )
        call rexcha(vother_qua    )
        call rexcha(w_vother_qua  )
        call rexcha(x0_qua        )
        call rexcha(y0_qua        )
        call rexcha(z0_qua        )
        call rexcha(massa_qua     )
        call rexcha(mezcla        )
        !
        ! Exchange of qua_reanut variables 
        !
        call iexcha(miinn_qua    )
        call iexcha(kfl_tiacc_qua)
        call iexcha(kfl_tisch_qua)

        call rexcha(cotol_qua)

        solve_sol => solve_qua(1:1)
        call soldef(1_ip)
        eigen_sol => eigen_qua(1:1)
        call eigdef(1_ip)
        !
        ! Exchange data read in qua_reaous
        !
        call iexcha(npp_inits_qua)           
        call iexcha(npp_iniso_qua)           
        do ji=1,nvarp_qua
           call iexcha(npp_stepi_qua(ji)) 
        end do
        call rexcha(pos_tinit_qua) 
        do ki=1,nvarp_qua
           do ji=1,nvart_qua
              call rexcha(pos_times_qua(ji,ki))
           end do
        end do
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
           if( ISLAVE .or. READ_AND_RUN() ) call Parall(2_ip)
        end if
     end do

     if( IMASTER .and. .not. READ_AND_RUN() ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','qua_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','qua_sendat',0_ip)     

     kfl_ptask = kfl_ptask_old

  case(2_ip)


     call Parall(27_ip)
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of qua_reabcs variables 
        !
        call iexcha(kfl_conbc_qua) 
        !call rexcha()    
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
           if( ISLAVE .or. READ_AND_RUN() ) call Parall(2_ip)
        end if
     end do

     if( IMASTER .and. .not. READ_AND_RUN() ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','qua_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','qua_sendat',0_ip)     

     if( ISLAVE ) call qua_membcs(1_ip)
     !
     ! Non-Constant boundary conditions 
     !
     if( kfl_conbc_qua == 0 .and. ISLAVE ) call qua_membcs(2_ip)

     if( INOTMASTER .or. .not. READ_AND_RUN() ) then
        !
        ! KFL_FIXNO_QUA
        !
        strin =  'KFL_FIXNO_QUA'
        call vocabu(NPOIN_INTE_1DIM,0_ip,0_ip)
        pari1 => kfl_fixno_qua(1,:)
        call Parall(300_ip)
        ! 
        ! BVESS_QUA
        !
        strre =  'BVESS_QUA'
        call vocabu(NPOIN_REAL_1DIM,0_ip,0_ip)
        parr1 => bvess_qua(1:,1)
        call Parall(300_ip)        
        !
        ! KFL_FIXBO_QUA
        !
        strin =  'KFL_FIXBO_QUA'
        call vocabu(NBOUN_INTE_1DIM,0_ip,0_ip)
        pari1 => kfl_fixbo_qua
        call Parall(300_ip)

        if( kfl_conbc_qua == 0 ) then
           !
           ! Non-constant bc
           !
           strin =  'BVESS_QUA'
           call vocabu(NPOIN_REAL_1DIM,0_ip,0_ip)
           parr1 => bvess_qua(1:,2)
           call Parall(300_ip) 

        end if
        !
        ! Master deallocates memory
        !
        if( IMASTER ) call qua_membcs(3_ip)

     end if

  case(3_ip)     
     !
     ! Exchange atomo
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     !call Parall(29_ip)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        nparc=nparc+50 
        if(parii==2.and.IMASTER) parch(1:50)          = atomo_qua(1) % wprob
        if(parii==2.and.ISLAVE)  atomo_qua(1) % wprob = parch(1:50)
        call iexcha( atomo_qua(1) % natoms )
        call iexcha( atomo_qua(1) % nespec )
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
           if( ISLAVE .or. READ_AND_RUN() ) call Parall(2_ip)
        end if
     end do

     if( IMASTER .and. .not. READ_AND_RUN() ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','qua_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','qua_sendat',0_ip)     
     !
     ! Allocatable arrays: ATOMO_QUA, ESPECIE_QUA
     !
     !call Parall(30_ip)
     if( ISLAVE ) then 
        call qua_memphy(1_ip)
        call qua_memphy(2_ip)
     end if

     kfl_ptask = kfl_ptask_old
     !
     ! Exchange especie
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        do ii = 1,atomo_qua(1) % natoms
           call iexcha( atomo_qua(1) % espe(ii)   )
           call iexcha( atomo_qua(1) % Tiespe(ii) )
           call iexcha( atomo_qua(1) % spin(ii)   )
           call iexcha( atomo_qua(1) % tipoPP(ii) )
           call rexcha( atomo_qua(1) % coorx(ii)  )
           call rexcha( atomo_qua(1) % coory(ii)  )
           call rexcha( atomo_qua(1) % coorz(ii)  )
        end do
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
           if( ISLAVE .or. READ_AND_RUN() ) call Parall(2_ip)
        end if
     end do

     if( IMASTER .and. .not. READ_AND_RUN() ) call Parall(2_ip)
print*,'b=',kfl_paral

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','qua_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','qua_sendat',0_ip)     
 

     kfl_ptask = kfl_ptask_old

     !
     ! Exchange especie
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     !call Parall(29_ip)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        do ii = 1,atomo_qua(1)%nespec
           nparc=nparc+100 
           if(parii==2.and.IMASTER) parch(1:50)              = especie_qua(ii) % wprob1
           if(parii==2.and.ISLAVE)  especie_qua(ii) % wprob1 = parch(1:50)
           if(parii==2.and.IMASTER) parch(51:100)            = especie_qua(ii) % wprob2
           if(parii==2.and.ISLAVE)  especie_qua(ii) % wprob2 = parch(51:100)
           call iexcha( especie_qua(ii) % atnumber )
           call iexcha( especie_qua(ii) % ncp      )
           call iexcha( especie_qua(ii) % ncl      )
           call iexcha( especie_qua(ii) % ncm      )
           call iexcha( especie_qua(ii) % nspin    )
           call iexcha( especie_qua(ii) % nlmax    )
           call iexcha( especie_qua(ii) % nlocal   )
           call rexcha( especie_qua(ii) % valencia )
           call iexcha( especie_qua(ii) % nrad     )
        end do

        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
           if( ISLAVE .or. READ_AND_RUN() ) call Parall(2_ip)
        end if
     end do

     if( IMASTER .and. .not. READ_AND_RUN() ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','qua_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','qua_sendat',0_ip)     
 

     kfl_ptask = kfl_ptask_old

     !
     ! Exchange especie
     !
     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     call Parall(29_ip)

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        do ii = 1,atomo_qua(1)%nespec
           do ji = 1,especie_qua(ii) % nrad
              call rexcha( especie_qua(ii) % radio(ji) )
              do ki = 1,especie_qua(ii) % nlmax
                 call rexcha( especie_qua(ii) % ppseu(ji,ki) )
                 call rexcha( especie_qua(ii) % ppphi(ji,ki)  )
              end do
           end do
           do ji = 1,100
              call rexcha( especie_qua(ii) % nocupa(ji)    )
           end do
        end do
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
           if( ISLAVE .or. READ_AND_RUN() ) call Parall(2_ip)
        end if
     end do

     if( IMASTER .and. .not. READ_AND_RUN() ) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','qua_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','qua_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','qua_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','qua_sendat',0_ip)     
 

     kfl_ptask = kfl_ptask_old

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine qua_sendat
