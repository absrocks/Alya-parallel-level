subroutine chm_updtsf(dtnew)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updtsf
  ! NAME 
  !    chm_updtsf
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
  real(rp),    intent(inout)  :: dtnew
  integer(ip)                 :: ispec,jspec,kspec,kcrit
  integer(ip)                 :: inisp,finsp,numsp
  real(rp)                    :: dt1,dt2,xfac1,dt1dt2,dt2dt1,s,a
  real(rp)                    :: dtacc,xfact,facc,dtfin,dtall,d
  real(rp)                    :: dtht,dttr,err1,err2
  integer(ip), save           :: ipass=0
  integer(4)                  :: istat
  integer(ip), save, pointer  :: lnerg_chm(:)
  real(rp),    save, pointer  :: energ_chm(:,:)
  real(rp),    save, pointer  :: xnerg_chm(:,:)
  real(rp)                    :: e3e1,e2e1
  real(rp)                    :: u1,u2,u3,dudt,d2udt2,dudtb,d2udt2b

  !----------------------------------------------------------------------
  !
  ! Allocate memory
  !
  !----------------------------------------------------------------------

  if( ipass == 0 ) then
     ipass = ipass + 1
     kcrit = 0
     kspec = 0
     allocate( lnerg_chm(nspec_chm),     stat=istat )
     allocate( energ_chm(nspec_chm+1,3), stat=istat )
     allocate( xnerg_chm(nspec_chm+1,3), stat=istat )
     do ispec = 1,nspec_chm+1
        energ_chm(ispec,1) = 0.0_rp
        energ_chm(ispec,2) = 0.0_rp
        energ_chm(ispec,3) = 0.0_rp
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Compute energy
  !
  !----------------------------------------------------------------------

  call chm_volume(energ_chm)

  if( ittim <= 4 ) then

     dtnew = dtold(1)

  else

     !----------------------------------------------------------------------
     !
     ! Target
     !
     !----------------------------------------------------------------------

     if( kfl_dttar_chm == -1 ) then
        inisp = 1
        finsp = 1
        xnerg_chm(1,1) = energ_chm(nspec_chm+1,1)
        xnerg_chm(1,2) = energ_chm(nspec_chm+1,2)
        xnerg_chm(1,3) = energ_chm(nspec_chm+1,3)
     else if( kfl_dttar_chm >= 1 ) then
        inisp = 1
        finsp = kfl_dttar_chm
        do ispec = 1,nspec_chm
           xnerg_chm(ispec,1) = energ_chm(ispec,1)
           xnerg_chm(ispec,2) = energ_chm(ispec,2)
           xnerg_chm(ispec,3) = energ_chm(ispec,3)
        end do
        call chm_sorten(3_ip,nspec_chm,xnerg_chm,lnerg_chm)
     end if

     !----------------------------------------------------------------------
     !
     ! Adaptive timestep
     ! Prediction of next timestep dtn+1 based on accuracy
     ! reached with previous timestep dtold.
     !     
     !----------------------------------------------------------------------

     dtacc  = 1.0e6_rp
     dt1    = dtold(1)  ! dt^n
     dt2    = dtold(2)  ! dt^n-1
     dt1dt2 = dt1/dt2   ! dt^n/dt^n-1
     dt2dt1 = dt2/dt1   ! dt^n-1/dt^n

     do ispec = inisp,finsp
        !
        ! Ratios
        !
        if( xnerg_chm(ispec,1) /= 0.0_rp ) then
           e3e1  = xnerg_chm(ispec,3) / xnerg_chm(ispec,1)
           e2e1  = xnerg_chm(ispec,2) / xnerg_chm(ispec,1)
           xfac1 = 1.0_rp - (1.0_rp+dt1dt2)*e2e1 + dt1dt2*e3e1 
        else
           e3e1  = 1.0_rp
           e2e1  = 1.0_rp
           xfac1 = 0.0_rp
        end if
        !
        ! Previous time step values
        !
        u1      = xnerg_chm(ispec,1)
        u2      = xnerg_chm(ispec,2)
        u3      = xnerg_chm(ispec,3)

        dudt    = (   u1*(dt2dt1+2.0_rp) -   u2*(dt1dt2+dt2dt1+2.0_rp) +   u3*dt1dt2) / (dt1+dt2)
        dudtb   = u1 * ( (dt2dt1+2.0_rp) - e2e1*(dt1dt2+dt2dt1+2.0_rp) + e3e1*dt1dt2) / (dt1+dt2)

        d2udt2  = (   u1*dt2dt1 -   u2*(1.0_rp+dt2dt1) +   u3) / (0.5_rp*dt2*(dt1+dt2))
        d2udt2b = u1 * ( 1.0_rp - (1.0_rp+dt1dt2)*e2e1 + dt1dt2*e3e1 ) / ( 0.5_rp * dt1 * ( dt1 + dt2 ) ) 

        !write(88,'(20(a,e16.8E3))') 'errors=',u1,' ',u2,' ',u3,' ',dudt,' ',dudtb,' ',d2udt2,' ',d2udt2b
        !write(89,'(20(a,f))')       'errors=',e2e1,' ',e3e1

        if( abs(xfac1) > zeror .and. abs(d2udt2b) > zeror ) then
           err1 = 0.5_rp * dt1 * ( dt1 + dt2 ) / xfac1        ! |e| / |d^2e/dt^2| (local truncation) 
           err2 =  ( dt2dt1 + (dt1dt2-dt2dt1)*e2e1 &          ! |de/dt| / |d^2e/dt^2| (high order terms)
                &  - dt1dt2*e3e1 ) / xfac1 * 0.5_rp * dt1
           dtht = sqrt(2.0_rp * epsht_chm * abs(err1)) 
           dttr = 2.0_rp * epstr_chm * abs(err2)        
        else
           dtht = strec_chm * dtold(1)
           dttr = strec_chm * dtold(1)
        end if
        !
        ! Criteria
        !
        if( kfl_dtcri_chm == 1 ) then
           dtall = dtht
           kcrit = 1
        else
           if( dtht < dttr ) then
              kcrit = 1
           else
              kcrit = 2
           end if
           dtall = min(dtht,dttr)
        end if
        if( dtall < dtacc ) then
           kspec = ispec
           dtacc = dtall
        end if

     end do
     !
     ! Calculation of accuracy reached with precious timestep dtold:
     ! facc = max(dt_trunc,dt_HT) / dtold
     !
     facc = dtacc / dtold(1)
     !
     ! Prediction of next time step dtn+1 to be used in the solver
     ! High damping: more stiff
     ! Low damping:  more smooth
     !
     a = facc
     d = dampi_chm
     if( a >= 1.0_rp ) then
        s = strec_chm
     else
        s = 1.0_rp/strec_chm
     end if
     xfact = (s-1.0_rp)*tanh( (a-1.0_rp)/(d*(s-1.0_rp)) ) + 1.0_rp

     dtnew = dtold(1) * xfact
     !
     ! Limit time step from below and above
     !
     if( dtnew < dtmin_chm ) then
        if( dtold(1) < dtmin_chm ) then 
           dtnew = dtold(1) * strec_chm
        else
           dtnew = dtmin_chm
        end if
     end if
     if( dtnew > dtmax_chm ) then
        if( dtold(1) > dtmax_chm ) then 
           dtnew = dtold(1) / strec_chm
        else
           dtnew = dtmax_chm
        end if
     end if
     !
     ! Check if final time step has been encompassed
     !
     dtfin = timef - cutim
     if( dtnew > dtfin ) dtnew = dtfin + zeror
 
     !----------------------------------------------------------------------
     !
     ! Write results in file
     !
     !----------------------------------------------------------------------

     if( INOTSLAVE ) then
        write(lun_times_chm,1) cutim,kcrit,kspec,dtnew,xfact,facc,dtmat_chm
        if( kfl_dttar_chm >= 1 ) &
             write(lun_time2_chm,2) (lnerg_chm(ispec),ispec=inisp,min(finsp,100_ip))
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Save previous values
  !
  !----------------------------------------------------------------------

  do ispec = 1,nspec_chm+1
     energ_chm(ispec,3) = energ_chm(ispec,2)
     energ_chm(ispec,2) = energ_chm(ispec,1)
  end do

1 format((1x,e12.6),2(1x,i7),10(1x,e12.6))
2 format(100(1x,i6))

end subroutine chm_updtsf
