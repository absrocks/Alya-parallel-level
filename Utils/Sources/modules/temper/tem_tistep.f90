subroutine tem_tistep
  !-----------------------------------------------------------------------
  !****f* Temper/tem_tistep
  ! NAME 
  !    tem_tittim
  ! DESCRIPTION
  !    This routine sets the time step
  ! USES
  ! USED BY
  !    tem_begite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod, only : kfl_adj_prob
  use def_temper
  use mod_ADR, only : ADR_time_strategy
  use mod_outfor, only : outfor
  implicit none

!!$  if(kfl_timco/=2) then
!!$     dtinv_tem=dtinv
!!$     if(kfl_stead_tem==1) dtinv_tem = 0.0_rp
!!$     if(kfl_timei_tem==0) dtinv_tem = 0.0_rp  
!!$    
!!$     if(kfl_tisch_tem==1) then
!!$        !
!!$        ! Trapezoidal rule: Euler iterations
!!$        !
!!$        pabdf_tem(1) = 1.0_rp
!!$        pabdf_tem(2) = -1.0_rp
!!$        if(ittim<=neule_tem) then
!!$           kfl_tiacc_tem=1
!!$        else
!!$           kfl_tiacc_tem=kfl_tiaor_tem
!!$        end if
!!$       
!!$        if(kfl_tiacc_tem==2) then
!!$           pabdf_tem(1) =  2.0_rp
!!$           pabdf_tem(2) = -2.0_rp
!!$           if (kfl_sgsac_tem==2) then
!!$              pabds_tem(1) =  2.0_rp
!!$              pabds_tem(2) = -2.0_rp
!!$           end if
!!$        end if
!!$        nbdfp_tem = 2
!!$        nbdfs_tem = 2
!!$     else if(kfl_tisch_tem==2) then
!!$        !
!!$        ! BDF scheme: increase integration order at each time step
!!$        !     
!!$        if (ittim<=neule_tem) then
!!$           kfl_tiacc_tem=1
!!$        else
!!$           kfl_tiacc_tem=min(kfl_tiaor_tem,ittim)
!!$        end if
!!$        call parbdf(kfl_tiacc_tem,pabdf_tem)        
!!$        call parbdf(min(kfl_sgsac_tem,kfl_tiacc_tem),pabds_tem)
!!$        nbdfp_tem = kfl_tiacc_tem+1
!!$        nbdfs_tem = min(nbdfp_tem, kfl_sgsac_tem+1)
!!$     end if
!!$  end if
  !
  ! Actualize time integration parameters: above is obsolete
  !
  if(kfl_adj_prob == 1 ) dtinv = 0.0_rp
  call ADR_time_strategy(ittim,dtinv,dtinv_old,ADR_tem)

  routp(1) = dtcri_tem
  routp(2) = 0.0_rp
  routp(3) = 0.0_rp
  ioutp(1) = kfl_timei_tem
  ioutp(2) = kfl_stead_tem
  call outfor(8_ip,lun_outpu,' ')

end subroutine tem_tistep
