subroutine chm_updtse(dtmin)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updtse
  ! NAME 
  !    chm_updtse
  ! DESCRIPTION
  !    This routine computes next timestep based on the accuracy reached
  !    with previous timestep
  ! USED BY
  !    chm_timste
  !***
  !----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  real(rp),    intent(inout)  :: dtmin
  integer(ip)                 :: ipoin,ispec,ipoi1,ipoi2,ii,kspec,kcrit
  real(rp)                    :: dt1,dt2,xfac1,xfac2,xfac3,xfac4,xfac5
  real(rp)                    :: dt1dt2,dt2dt1,s,a
  real(rp)                    :: dta,dtb,e1,e2,e3,dtacc,xfact,dtac1,facc
  real(rp)                    :: dtht,dttr,e3e1,e2e1,err1,err2,fx,gx,hx,vx
  integer(ip), save           :: ipass=0
  integer(4)                  :: istat
  real(rp),    save           :: energ_chm(3)

  if( ipass == 0 ) then
     ipass        = ipass + 1
     energ_chm(1) = 0.0_rp
     energ_chm(2) = 0.0_rp
     energ_chm(3) = 0.0_rp
  end if
  !
  ! Compute energy
  !
  energ_chm(1) = 0.0_rp
  dtht         = 1.0e6_rp
  dttr         = 1.0e6_rp

  if( INOTMASTER ) then     
     do ispec = 1,nspec_chm        
        do ipoin = 1,npoin
           energ_chm(1) = energ_chm(1) + conce(ipoin,ispec,3) * conce(ipoin,ispec,3)
        end do
     end do
     energ_chm(1) = energ_chm(1) / real(npoin)
  end if
  energ_chm(1) = postp(1) % veset(1,1)


  call pararr('SUM',0_ip,1_ip,energ_chm)

  if( ittim <= 4 ) then

     dtmin = dtold(1)

  else
     !
     ! Adaptive timestep
     ! Prediction of next timestep dtn+1 based on accuracy
     ! reached with previous timestep dtold.
     !     
     if( INOTMASTER ) then

        dta   = dtmin
        dtb   = 100.0_rp * dtmin
        dtacc = 1.0e6_rp
        dt1   = dtold(1)
        dt2   = dtold(2)
        xfac1 = dt1*dt2*(dt1+dt2)
        xfac2 = dt2*(dt2+2.0_rp*dt1)
        xfac3 = (dt1+dt2)*(dt1+dt2)
        xfac4 = dt1*dt1
        xfac5 = dt1+dt2
        dt1dt2 = dt1/dt2
        dt2dt1 = dt2/dt1
        !
        ! Ratios
        !
        if( energ_chm(1) /= 0.0_rp ) then
           e3e1  = energ_chm(3) / energ_chm(1)
           e2e1  = energ_chm(2) / energ_chm(1)
        else
           e3e1 = 1.0_rp
           e2e1 = 1.0_rp
        end if
        !
        ! Previous time step values
        !
        xfac3 = ( dt2dt1 + e3e1 ) - e2e1 * ( 1.0_rp + dt2dt1 ) 
        if( abs(xfac3) > zeror ) then
           err1 = 0.5_rp * dt2 * ( dt1 + dt2 ) / xfac3
           err2 = ( dt2dt1 + 2.0_rp + dt1dt2 * e3e1 - e2e1 * ( dt1dt2 + dt2dt1 + 2.0_rp ) ) &
                & / xfac3 * 0.5_rp * dt2
           dtht = sqrt(2.0_rp * epsht_chm * abs(err1))  !  |e| / |d^2e/dt^2|
           !dttr(ispec) = 2.0_rp * epstr_chm * abs(err2)        !  |de/dt| / |d^2e/dt^2|
           !dttr(ispec) = 1.0e6_rp
        end if
        !
        ! Look for maximum over subdomains
        !
        call pararr('MIN',0_ip,1_ip,dtacc)
        !
        ! Calculation of accuracy reached with precious timestep dtold:
        ! facc=max(dt_trunc,dt_HT)/dtold
        !
        dtacc = min( dtht , dttr )
        facc  = dtacc / dtold(1)
        !
        ! Prediction of next time step dtn+1 to be used in the solver:
        ! dtn+1=dtold * min(strech,ln(1+(damp-1)*facc)/ln(damp))
        ! 
        ! strech: strech factor to limit timestep increase such thate dtnew <= P2.dtold.
        ! strech_m to be defined by the user in INPUT. strech >1. (eg: strech=1.05 means 5%).
        ! damp: damping factor to avoid oscillations of the timestep.
        ! Given by user in INPUT. damp > 1.0.
        !
        !if( facc > 0.5_rp ) then
        if( 1 == 0 ) then
           xfact = min(strec_chm,log(1.0_rp + (dampi_chm-1.0_rp)*facc) / log(dampi_chm))
           xfact = max(xfact,1.0_rp/strec_chm)
        else
           a     = strec_chm
           s     = ( a**2+1.0_rp)/(2.0_rp*a)
           xfact = 0.5_rp*(strec_chm-1.0_rp/strec_chm)*tanh(1.0_rp*(facc-s))+s
        end if
        dtmin = dtold(1) * xfact
        dtmin = max(1.0e-12_rp,dtmin)

     end if
     if( INOTSLAVE ) then
        write(lun_times_chm,1) ittim,dtmin,xfact,&
             log(1.0_rp+(dampi_chm-1.0_rp)*facc)/log(dampi_chm),facc,energ_chm(1)
     end if

  end if

  energ_chm(3) = energ_chm(2)
  energ_chm(2) = energ_chm(1)

1 format(1(1x,i7),10(1x,e12.6))

end subroutine chm_updtse
