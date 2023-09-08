subroutine got_sendat(order)
!-----------------------------------------------------------------------
!****f* gotita/got_sendat
! NAME
!    got_sendat
! DESCRIPTION
!    This routine exchange input data 
! USES
! USED BY
!    got_turnon
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_domain
  use def_gotita
  use def_inpout
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,ki,dummi
  integer(4)              :: istat

  select case (order)

  case(1)     
     !
     ! Exchange data read in got_reaphy, got_reanut and got_reaous
     !
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        !
        ! Exchange of got_reaphy variables 
        !
        call iexcha(kfl_diffu_got)
        call iexcha(kfl_difun_got)
        call iexcha(kfl_forme_got)
        call iexcha(kfl_probl_got)
        call iexcha(kfl_timei_got)
        call iexcha(kfl_timec_got)
        call iexcha(kfl_timem_got)
        call iexcha(kfl_velfu_got)   
        call rexcha(ddrop_got)
        call rexcha(deair_got)
        call rexcha(densi_got)
        do ji=1,3
           call rexcha(diffu_got(ji))
        end do
        do ji=1,3
           call rexcha(gravi_got(ji))
        end do
        call rexcha(grnor_got)
        call rexcha(leinf_got)
        call rexcha(muair_got)
        call rexcha(veair_got)
        !
        ! Exchange of got_reanut variables 
        !
        call iexcha(kfl_algor_got)
        call iexcha(kfl_artif_got)
        call iexcha(kfl_coupl_got)
        call iexcha(kfl_dttyp_got)
        call iexcha(kfl_ellen_got)
        call iexcha(kfl_linea_got)
        call iexcha(kfl_normc_got)
        call iexcha(kfl_penal_got)
        call iexcha(kfl_sgsco_got)
        call iexcha(kfl_sgsti_got)
        call iexcha(kfl_shocc_got)
        call iexcha(kfl_shocm_got)
        call iexcha(kfl_staty_got)
        call iexcha(kfl_taust_got)
        call iexcha(kfl_tiacc_got)
        call iexcha(kfl_weigc_got)
        call iexcha(kfl_weigm_got)
        call iexcha(itart_got)
        call iexcha(itshc_got)
        call iexcha(itshm_got)
        call iexcha(mibgs_got)
        call iexcha(miinn_got)
        call iexcha(misgs_got)
        call iexcha(neule_got)
        call iexcha(npica_got)
        call rexcha(artif_got(1))
        call rexcha(artif_got(2))
        call rexcha(cotol_got)
        call rexcha(cutof_got)
        call rexcha(penal_got)
        call rexcha(relax_got)
        call rexcha(relgs_got)
        call rexcha(relsg_got)
        call rexcha(safet_got)
        call rexcha(shock_got)
        call rexcha(sstol_got)
        call rexcha(staco_got(1))
        call rexcha(staco_got(2))
        call rexcha(staco_got(3))
        call rexcha(tobgs_got)
        call rexcha(tosgs_got) 
        do jr=1,3
           call rexcha(xmaxi_got(jr))
        end do
        do jr=1,3
           call rexcha(xmini_got(jr))
        end do
        solve_sol => solve(1:)
        call soldef(1_ip)
        solve_sol => solve(2:)
        call soldef(1_ip)
        solve_sol => solve(3:)
        call soldef(1_ip)
        !
        ! Exchange data read in got_reaous
        !
        call posdef(1_ip,dummi)
        call iexcha(kfl_exacs_got)
        call iexcha(npp_bound_got)
        call rexcha(pos_cutof_got)
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','got_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','got_sendat',parre)
           if(kfl_paral>=1) call Parall(two)
        end if
     end do

     if(kfl_paral==0) call Parall(two)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','got_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','got_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','got_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','got_sendat',0_ip)

     if(kfl_probl_got/=3.and.kfl_velfu_got==-1) then
        !
        ! Exchange air velocity
        !        
        if(kfl_paral/=0) call got_memphy(1_ip)
        call vocabu(NPOIN_REAL_2DIM,ndime,0_ip)
        parr2 => veloc_got(:,:) 
        call Parall(300_ip)

     else if(kfl_probl_got==3) then
        !
        ! Exchange droplet velocity
        !
        if(kfl_paral/=0) call got_memphy(2_ip)
        call vocabu(NPOIN_REAL_2DIM,ndime,0_ip)
        parr2 => vdrop(:,:,1) 
        call Parall(300_ip)

     end if

  case (2)
       !
       ! Exchange data read in got_reabcs
       !
       if(kfl_paral>=1) call got_membcs(one)
       !
       ! KFL_FIXNO_GOT
       !
       call vocabu(NPOIN_INTE_1DIM,0_ip,0_ip)
       pari1 => kfl_fixno_got
       call Parall(300_ip)
       !
       ! BVESS_GOT
       !
       call vocabu(NPOIN_REAL_2DIM,ndofn_got(3),0_ip)
       parr2 => bvess_got(:,:) 
       call Parall(300_ip)
       !
       ! Deallocate master's memory
       !
       if(kfl_paral==0) call got_membcs(two)

  end select

end subroutine got_sendat
