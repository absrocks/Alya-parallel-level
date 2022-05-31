!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_cvgunk.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Convergence
!> @details Convergence
!> @} 
!-----------------------------------------------------------------------

subroutine neu_cvgunk(itask)

  use def_parame
  use def_master
  use def_domain
  use def_neutro
  use mod_outfor, only : outfor
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: idirection,ienergy
  integer(ip), save       :: ipass=0
  real(rp),    save       :: cpuit_neu=0.0_rp
  real(rp)                :: time1,ritl2
  real(rp),    save       :: residual_neu

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  select case( itask )

  case ( 1_ip )
     !
     ! L2 residual
     !
     if( current_energy_neu == 1 .and. current_direction_neu == 1 ) residual_neu = 0.0_rp
     call residu(2_ip,one,one,unkno,neutr(current_energy_neu,current_direction_neu,1:npoin,1),one,one,one,1.0_rp,ritl2)     
     residual_neu = residual_neu + ritl2**2 
     
  case ( 4_ip )

     residual_neu = sqrt(residual_neu)
     if( residual_neu < cotol_neu .or. itinn(modul) >= miinn_neu ) kfl_goite_neu = 0
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        call cputim(time1)
        if( ipass == 0 ) then
           time1=time1 - cpu_initi
        else
           time1=time1 - cpuit_neu
        end if
        if( ipass == 0 .and. kfl_rstar /= 2 ) write(momod(modul)%lun_conve,100)
        write(momod(modul)%lun_conve,101) ittim,itcou,itinn(modul),cutim,residual_neu
        call cputim(cpuit_neu)
        flush(momod(modul)%lun_conve)
        ipass = 1
     end if

  case ( 2_ip )
     !
     ! Check convergence of the outer iterations:
     ! || T(n,i,*) - T(n,i-1,*)|| / ||T(n,i,*)||
     !    
     do idirection = 1,num_directions_neu
        do ienergy = 1,num_energies_neu
           call residu(2_ip,one,one,neutr(ienergy,idirection,1:npoin,1),neutr(ienergy,idirection,1:npoin,2),one,one,one,1.0_rp,ritl2)
           resid_neu = resid_neu + ritl2**2
        end do
     end do
     resid_neu = sqrt(resid_neu)

  case ( 3_ip )
     !
     ! Check residual of the time iterations:
     ! || T(n,*,*) - T(n-1,*,*)|| / ||T(n,*,*)||
     !
!!$     call residu(2_ip,one,one,therm(1,1),therm(1,min(3_ip,ncomp_neu)),one,one,one,1.0_rp,rineu)
!!$     if(rineu<=sstol_neu) then
!!$        kfl_stead_neu = 1
!!$        call outfor(28_ip,momod(modul)%lun_outpu,' ')
!!$     end if
!!$     !
!!$     ! Low-Mach model
!!$     !
!!$     if( INOTSLAVE .and. kfl_regim_neu>=3 ) then
!!$        if( jpass == 0 ) then
!!$           jpass = 1
!!$           write(lun_lmach_neu,400)
!!$        end if
!!$        write(lun_lmach_neu,'(2x,10(2x,e12.6))') cutim,prthe(1),prthe(1)/prthe(4),dpthe,xmass_neu
!!$     end if

  end select

  !
  ! Formats
  !
100 format('# --| ALYA Convergence '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       &   '# --| 4. Current time      5. Radiation          ' ,//,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5')
101 format(4x,i9,2x,i9,2x,i9,11(2x,e12.6))

end subroutine neu_cvgunk

