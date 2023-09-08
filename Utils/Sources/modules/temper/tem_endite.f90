subroutine tem_endite(itask)
!-----------------------------------------------------------------------
!****f* Temper/tem_endite
! NAME 
!    tem_endite
! DESCRIPTION
!    This routine checks convergence and performs updates of the
!    temperature  at:
!    - itask=1 The end of an internal iteration
!    - itask=2 The end of the internal loop iteration
! USES
!    tem_cvgunk
!    tem_updunk
! USED BY
!    tem_doiter
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain,    only : elmar,npoin
  use def_temper
  use mod_ker_proper
  use mod_ADR,       only : ADR_manufactured_error
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: itask,ipoin,ivalu
  real(rp)    :: dummr,tenew,cploc(6,2)

  select case(itask)

  case(1)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || T(n,i,j) - T(n,i,j-1)|| / ||T(n,i,j)||) and update unknowns:
     !  T(n,i,j-1) <-- T(n,i,j) 
     !
     call tem_updunk(7_ip)    ! Compressible flow update rho
     call tem_updunk(10_ip)   ! Cut off undershoots
     call tem_cvgunk(1_ip)    ! Residual:   ||UNKNO(:)-TEMPE(:,1)||
     call tem_updunk(3_ip)    ! Relaxation and update:     TEMPE(:,1)=UNKNO

     if(kfl_regim_tem>=3) then
        ! 
        ! If low Mach Updates thermodynamic pressure 
        !
        call tem_updthe(1)
     end if
     !
     ! solves subgrid scales
     !  
     call tem_solsgs()
     ! 
     ! If low Mach with sgs updates thermodynamic pressure again
     !
     if(kfl_regim_tem>=3.and. kfl_sgsti_tem /= 0) then
        call tem_updthe(1)
     end if
 
  case(2)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || T(n,i,*) - T(n,i-1,*)|| / ||T(n,i,*)||) and update unknowns:
     !  T(n,i-1,*) <-- T(n,i,*) 
     !
     call livinf(16_ip,' ',itinn(modul))
     call tem_cvgunk(2_ip)            ! Residual: ||TEMPE(:,2)-TEMPE(:,1)||
     call tem_updunk(4_ip)            ! Update:   TEMPE(:,2) = TEMPE(:,1)

     if (kfl_regim_tem == 4) then
       ! 
       ! Re-compute the enthalpy with the new values of C_p to ensure the correct T at the boundary 
       !
       if (kfl_regim_tem==4 .and. kfl_plepp_tem /= 4) call tem_calcEnthalpyBC() 
       ! 
       ! Compute current temperature from enthalpy 
       !
       call tem_comtem
       ! 
       ! Store current temperature (:,2) <= (:,1)
       ! 
       call tem_store_tempe(one) 

     end if

     !!call tem_updhfl()     ! Update:   TFLUX= Low-Mach heat flux contribution

     ! ********************
     !   Test for coupling 
     ! call tem_temexch(2_ip)
     ! call tem_temexch(1_ip)
     ! ********************
     if( kfl_exacs_tem /= 0 .and. kfl_timei_tem == 0 ) then
        call ADR_manufactured_error(ADR_tem,ittim,cutim,therm)
     end if
     
     if (kfl_cos_opt == 1) call tem_costcal()
     if (kfl_adj_prob == 1) call tem_senscal()

  end select

end subroutine tem_endite
