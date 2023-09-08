!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_cvgunk.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Performs several convergence checks for the porous equation
!> @details Performs several convergence checks for the porous equation
!> @} 
!------------------------------------------------------------------------
subroutine por_cvgunk(itask)
  use      def_parame
  use      def_master
  use      def_domain
  use      def_porous
  use mod_outfor, only : outfor
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0
  real(rp),    save       :: cpuit_por=0.0_rp
  real(rp)                :: ripor(2),time1,ritl2
  real(rp)                :: prmin,prmax  ! in temperature they are part of def temper but I see no sense if they are only used here
  real(rp)                :: samin,samax

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  select case(itask)

  case(1)
     !
     ! Check convergence of the inner iterations:
     ! || P(n,i,j) - P(n,i,j-1)|| / ||P(n,i,j)||
     !
     if(kfl_normc_por == 3) then  
        ripor=solve(1)%resin   ! Algebraic residual 
     else  ! see how we should deal with this - should we check only convergence for the pressure or also for the water sat
        call residu(kfl_normc_por,one,one,unkno,press,one,one,one,1.0_rp,ritl2)     
        ripor(1)=ritl2                ! L2 residual
     end if
     
     if( (ripor(1)<cotol_por).or.(itinn(modul)>=miinn_por)) kfl_goite_por = 0  ! the pressure is in charge of stopping the iteration
     !
     ! Compute min and max of the porous unknowns  -- see how to treat this
     !
     call minmax(one,npoin,zero,unkno,prmin,prmax)
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        call cputim(time1)
        if(ipass==0) then
           time1 = time1-cpu_initi
        else
           time1 = time1-cpuit_por
        end if
        if(ipass==0.and.kfl_rstar/=2) write(momod(modul)%lun_conve,100)
        write(momod(modul)%lun_conve,101) ittim,itcou,itinn(modul),cutim,ripor,resgs_por(1),&
             prmin,prmax,time1
        call cputim(cpuit_por)
        flush(momod(modul)%lun_conve)
        ipass=1
     end if

  case(2)
     !
     ! Check convergence of the outer iterations:   ! for the moment only checkin press -- satur???
     ! || P(n,i,*) - P(n,i-1,*)|| / ||P(n,i,*)||
     !
     call residu(2_ip,one,one,press(1,1),press(1,2),one,one,one,1.0_rp,resid_por(1))
     call residu(2_ip,one,one,wasat(1,1),wasat(1,2),one,one,one,1.0_rp,resid_por(2))

  case(3)
     !
     ! Check residual of the time iterations:     ! for the moment only checkin press -- satur???
     ! || P(n,*,*) - P(n-1,*,*)|| / ||P(n,*,*)||
     !
     call residu(2_ip,one,one,press(1,1),press(1,3),one,one,one,1.0_rp,ripor(1))
     call residu(2_ip,one,one,press(1,1),wasat(1,3),one,one,one,1.0_rp,ripor(2))
     if ( (ripor(1)<=sstol_por) .and. (ripor(2)<=sstol_por) ) then
        kfl_stead_por = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ')
     end if

  case(4)
     !
     ! Check convergence of the inner iterations:
     ! || S(n,i,j) - S(n,i,j-1)|| / ||S(n,i,j)||
     !
     if(kfl_normc_por==3) then      ! While we are doing saturation explicitly this does not make sense - so I use L2 norm
!        ripor=solve(1)%resin   ! Algebraic residual 
        call residu(2_ip,one,one,unkno,wasat,one,one,one,1.0_rp,ritl2)     
        ripor(2) = ritl2                ! L2 residual
     else  
        call residu(kfl_normc_por,one,one,unkno,wasat,one,one,one,1.0_rp,ritl2)     
        ripor(2) = ritl2                ! L2 residual
     end if
     
!     if((ripor(2)<cotol_por).or.(itinn(modul)>=miinn_por)) kfl_goite_por = 0 ! the pressure is in charge of stopping the iteration
     !
     ! Compute min and max of the Saturation
     !
     call minmax(one,npoin,zero,unkno,samin,samax)
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        call cputim(time1)
        if(ipass==0) then
           time1=time1-cpu_initi
        else
           time1=time1-cpuit_por
        end if
        if(ipass==0.and.kfl_rstar/=2) write(momod(modul)%lun_conve,100)
        write(momod(modul)%lun_conve,101) ittim,itcou,itinn(modul),cutim,ripor(1),ripor(2),&
             samin,samax,time1
        call cputim(cpuit_por)
        flush(momod(modul)%lun_conve)
        ipass=1
     end if

  end select
  !
  ! Formats
  !
100 format('# --| ALYA Convergence '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       &   '# --| 4. Current time      5. Pressure           6. Satuarion         ' ,/,& 
       &   '# --| 7. Min. pressure     8. Max. pressure      9. Min. Satuartion  ' ,//,&
       &   '# --| 10. Max. Saturation  11. Elapsed CPU Time  ' ,//,&
       &   '$ ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9','            10','            11')
101 format(4x,i9,2x,i9,2x,i9,13(2x,e13.6))

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine por_cvgunk

