subroutine Begste
  !-----------------------------------------------------------------------
  !****f* doiter
  ! NAME 
  !    doiter
  ! DESCRIPTION
  !    This routine sets some variable before solving a time step
  !    equations.
  ! USES
  !
  !
  ! USED BY
  ! Wind
  !***
  !-----------------------------------------------------------------------

  use      def_master
  implicit none
  real(8)      :: htime,limit, perio

  istep  = istep +1
  limit =1.0e-4
  if (istep.eq.1) safet = dtinv
  if (kfl_local) then
     dtinv = safet
     limit = 5.0e-2
  end if
  if (kfl_canop) limit = limit*3.0
  if (istep.gt.1.and..not.kfl_trtem) then
     !              dtinv=max(dtinv*0.7d0,2.0e-4)
     dtinv= max (dtinv*0.9d0,limit )
     if (kfl_local)   safet = dtinv 
  end if

  ctime = ctime +1.0d0/dtinv
  if (istep.eq.1.or.kfl_hflx_wall) tewal = tempe(1,1) 
  if (kfl_trtem) then      ! Transient temperature
     htime =ctime/3600.0d0 ! Current time in hours
     !
     ! Wall temperature 'tewal' evolution calculation
     !
     if (kfl_case.eq.1) then        ! GABLS2
        call gabls2_begste(htime)
        !
     else if (kfl_case.eq.2) then   ! Diurnal cycle
        call cycle_begste(htime)
        !
     else if (kfl_case.eq.3) then   ! GABLS3 (adds mesoscale tendencies)
        call gabls3_begste
        !
     else if (kfl_case.eq.4) then   ! GABLS1  ! May is now by defauls
        tewal =  265.0 - htime*0.25d0
        if (istep.eq.1)     teref = tewal
     else if (istep==1) then ! no prescribed case! 
       teref = tewal
  !        call runend('Begste: Transient THERMAL case is incorrect or not implemented')

     else ! no prescribed case! example: radiation case 
!        call interp_wallte(tewal)  ! only for prescribed wall temperature 
!        tewal = tewal +265.0 - 285.15
!        teref = tewal
!!$       perio = 24.0*3600.0
!!$       tewal = 265.0 - 8.0*sin(2.0*pi/perio*ctime)
!        print *, 'tewal',ctime/3600.0, tewal
       end if
  end if
end subroutine begste
