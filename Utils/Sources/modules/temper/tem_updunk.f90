!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_updunk.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   This routine performs several types of updates for the thermal variable (T or h)
!> @details Solution updates:
!>        PREPROCESSING
!>          tem_iniunk (itask=8)    ..................... if Diffusion problem
!>          tem_iniunk (itask=9)    ..................... if Diffusion problem
!>          tem_iniunk (itask=6)    ..................... if Restart
!>
!>        TIME LOOP
!>          do time
!>             tem_begste (itask=1) ..................... (:,2) <= (:,3)
!>             do outer
!>                tem_begite (itask=2) .................. (:,1) <= (:,2)
!>                do inner
!>                   tem_endite (itask=7 , inner loop) .. Update density (compressible flow)
!>                   tem_endite (itask=10, inner loop) .. Cut off undershoots
!>                   tem_endite (itask=3 , inner loop) .. (:,1) <= UNKNO
!>                end do
!>                tem_endite (itask=4, outer loop) ...... (:,2) <= (:,1)
!>             end do
!>             tem_endste (itask=5) ..................... (:,3) <= (:,1)
!>          end do
!> @}
!-----------------------------------------------------------------------

subroutine tem_updunk(itask)

  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod, only : gasco
  use mod_ADR,    only : ADR_end_time_step
  use mod_ADR,    only : ADR_begin_inner_iteration
  use mod_ADR,    only : ADR_end_inner_iteration
  use mod_ADR,    only : ADR_begin_time_step
  use mod_ADR,    only : ADR_after_restart
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ielem,igaus,itime,pgaus,ipoin,icomp
  real(rp)                :: rela1_tem

  if( INOTMASTER ) then

     select case (itask)

     case(1_ip)
        !
        ! (:,2) <= (:,3): Initial guess for outer iterations
        !
        icomp = min(3_ip,ncomp_tem)
        call ADR_begin_time_step(ADR_tem,therm) 
        if( kfl_regim_tem == 4 ) then
           do ipoin=1,nunkn_tem
              tempe(ipoin,2) = tempe(ipoin,icomp)
           end do
        endif

     case(2_ip)
        !
        ! (:,1) <= (:,2): Initial guess for inner iterations
        !
        call ADR_begin_inner_iteration(ADR_tem,therm) 

        do ipoin=1,nunkn_tem
           unkno(ipoin) = therm(ipoin,1)
        end do
        if( kfl_regim_tem == 4 ) then
           do ipoin = 1,nunkn_tem
              tempe(ipoin,1) = tempe(ipoin,2)
           end do
        end if

     case(3_ip)
        !
        ! Assign T(n,i,j-1) <-- T(n,i,j), update of the temperature
        !
        if(postp(1)%npp_stepi(8)>0) then
           do ipoin=1,nunkn_tem
              teold_tem(ipoin) = therm(ipoin,1)
           end do
        end if
        if(relax_tem==1.0_rp) then
           do ipoin=1,nunkn_tem
              therm(ipoin,1) = unkno(ipoin)
           end do
        else
           rela1_tem=1.0_rp-relax_tem
           do ipoin=1,nunkn_tem
              therm(ipoin,1) = relax_tem*unkno(ipoin)+rela1_tem*therm(ipoin,1)
           end do
        end if

     case(4_ip)
        !
        ! (:,2) <= (:,1): End of inner iteration
        !        
        call ADR_end_inner_iteration(ADR_tem,therm) 

     case(5_ip)
        !
        ! (:,3) <= (:,1): End of time step
        !        
        call ADR_end_time_step(ADR_tem,therm) 

     case(6_ip) 
        !
        ! (:,1) <= (:,3)
        ! 
        call ADR_after_restart(ADR_tem,therm) 

     case(7_ip)
        !
        ! Update density
        !
        if(kfl_regim_tem==1) then
           !
           ! Compressilbe flow: rho=P/(RT)
           !
           do ipoin=1,nunkn_tem
              densi(ipoin,1) = press(ipoin,1)/(gasco*tempe(ipoin,1))
           end do
        end if

     case(8_ip)
        !
        ! Solve initial problem
        !
        do ipoin=1,nunkn_tem
           unkno(ipoin) = therm(ipoin,1)
        end do

     case(9_ip)
        !
        ! Solve initial problem
        !
        icomp = min(3_ip,ncomp_tem)
        do ipoin=1,nunkn_tem
           therm(ipoin,1)     = unkno(ipoin) 
           therm(ipoin,icomp) = unkno(ipoin) 
        end do

     case(10_ip)
        !
        ! Prevent undershoots
        !
        if (kfl_negat_tem==1) then
           do ipoin = 1,nunkn_tem
              if( unkno(ipoin) < negat_tem ) then
                 unkno(ipoin) = negat_tem
              end if
           end do
        endif
        if (kfl_posit_tem==1) then
           do ipoin = 1,nunkn_tem
              if( unkno(ipoin) > posit_tem ) then
                 unkno(ipoin) = posit_tem
              end if
           end do
        endif

     end select

  end if

end subroutine tem_updunk
