subroutine readinp
  !********************************************************
  !*
  !*    Reads the properties from the input file
  !*
  !********************************************************
  use KindType
  use InpOut
  use def_Master
  implicit none
  !
  integer(ip)           :: istat,ivoid(1)
  real   (rp)           :: rvoid(mcuts)
  character(len=s_mess) :: message,cvoid

  !
  !*** initializations
  !
  dtinv           = 1.0d0     !Time step in seconds (then it is transformed to 1/dt)
  kfl_thcou       = .false.   !Thermal coupling
  kfl_trtem       = .false.   !Transient temperature
  kfl_canop       = .false.   !Canopi model on/off
  lmoni           = 0.0_rp    !Monin-Obukhov length (L)
  ztsbl           = 0.0_rp
  tewal           = teref
  ncuts           = 0
  !Models
  kfl_thmod       = 0         !No thermal coupling model  (default)
  kfl_canmo       = 0         !Sogachev's canopy model (default)
  kfl_candi       = 0         !uniform canopy distribution (default)
  kfl_abl2        = .false.   !ABL wall law
  kfl_logva       = .false.   !Logarithmic variables
  kfl_local       = .false.   !Local time step
  !Mesoscale tendencies
  kfl_thadv       = .false.   !Potential temperature advection
  kfl_momadv      = .false.   !Momentum components advection
  kfl_pressgr     = .false.   !Pressure gradient from mesoscale 
  !Nudging to mesoscale data
  kfl_mom_nudging = .false.   !for momentum
  kfl_tem_nudging = .false.   !for potential temperature
  !netCDF 
  kfl_netCDF      = .false.   !Activate netCDF output
  append_netCDF   = .false.   !Append outputs to existing netCDF file (for restart option)
  !Other
  kfl_topco_vel   = 0         !Top wind velocity free=0 or prescribed=1
  kfl_case        = 0         !Implemented thermally coupled cases
  !Restart
  restart_in      = .false.
  restart_out     = .false.
  heica =0.0d0
  damping         =0.0_rp
  z_damping       =0.0001_rp
  !
  !*** Writes to the log file
  !
  write(*,1)
  if(out_screen) write(*,1)
1 format(/,'---> Reading input data file...',/)

  !
  !*** Reads the GENERAL block
  ! 
  call get_input_int(finp,'GENERAL','model',ivoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  kfl_model = ivoid(1)

  call get_input_cha(finp,'GENERAL','logva',cvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  !  if(istat.lt.0) call runend(message)
  if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.&
       TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') then
     kfl_logva = .true. 
     print *, 'LOG VARIABLES ON!!'
  end if
  call get_input_cha(finp,'GENERAL','ABL2 ',cvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  !  if(istat.lt.0) call runend(message)
  if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.&
       TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') then
     kfl_abl2 = .true.
     print *, 'ABL2 wall model!!'
  end if

  call get_input_cha(finp,'GENERAL','LOCAL',cvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  !  if(istat.lt.0) call runend(message)
  if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.&
       TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') then
     kfl_local = .true.
     print *, 'Local time step!!'
  end if
  !
  call get_input_int(finp,'GENERAL','nstep',ivoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  nstep = ivoid(1)
  !
  call get_input_int(finp,'GENERAL','miite',ivoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  miite = ivoid(1) 
  !
  dtinv = 2.0d0 ! default time step
  call get_input_rea(finp,'GENERAL','timst',rvoid,1,istat,message)  
  if(istat.ne.0) call wriwar(message) 
  if (istat.eq.0)   dtinv = rvoid(1)
  dtinv = 1.0d0/dtinv

  !
  !*** Reads the RESTART block
  !
  call get_input_cha(finp,'RESTART','write',cvoid,1,istat,message)
  if(istat.lt.0) then
     call wriwar(message)
  elseif (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     restart_out = .true.
     write(*,*) 'RESTART: YES write'
  end if
  !
  if (restart_out) then
     call get_input_int(finp,'RESTART','freq',ivoid,1,istat,message)
     if(istat.gt.0) then
        call wriwar(message)
     elseif(istat.lt.0) then
        call runend(message)
     else
        freq_rst = ivoid(1)
!!$        if(mod(freq_rst,1.0d0/dtinv)/=0) then
!!$           call runend('netCDF output frequency is not a multiple of time step')
!!$        end if
     end if
  end if
  !
  call get_input_cha(finp,'RESTART','read',cvoid,1,istat,message)
  if(istat.lt.0) then
     call wriwar(message)
  elseif (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
     restart_in = .true.
     write(*,*) 'RESTART: YES reads Initial Condition'
  end if
  !
  !
  !*** Reads the PHYSICAL block
  !
  call get_input_rea(finp,'PHYSICAL','fcori',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  fcori = rvoid(1)

  call get_input_rea(finp,'PHYSICAL','vegeo',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  vegeo = rvoid(1)

  call get_input_cha(finp,'PHYSICAL','vegeo_file',cvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
!  if(istat.lt.0) call runend(message)
  vegeo_file = TRIM(cvoid)
  print *, 'vegeo_file =', TRIM(vegeo_file)

  
  call get_input_rea(finp,'PHYSICAL','ustar',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  ustar = rvoid(1) 
  call get_input_rea(finp,'PHYSICAL','l_max',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  l_max = rvoid(1)
  if (abs(fcori).lt.1.0d-9) then  
     l_max=1.0d50
  elseif (abs(l_max).lt.0.01) then
     l_max = 0.00027*vegeo/abs(fcori)
  end if
  print *, 'l_max=', l_max

  call get_input_rea(finp,'PHYSICAL','rough',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  rough = rvoid(1)
  !
  !*** Reads the MESH block
  !
  call get_input_rea(finp,'MESH','ztop',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  length = rvoid(1)

  call get_input_rea(finp,'MESH','hele1',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  dz1 = rvoid(1)

  call get_input_int(finp,'MESH','nelem',ivoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  nelem = ivoid(1)

  !
  !*** Reads the BOUNDARY block
  !
  call get_input_int(finp,'BOUNDARY','bouco',ivoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  kfl_bouco_vel = ivoid(1)

  call get_input_int(finp,'BOUNDARY','topco',ivoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  ! if(istat.lt.0) call runend(message)
  kfl_topco_vel = ivoid(1)
  !
  !*** Reads the NUMERICAL block
  !
  call get_input_rea(finp,'NUMERICAL','dwall',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  dwall = rvoid(1)

  call get_input_rea(finp,'NUMERICAL','toler',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  toler(1) = rvoid(1)
  toler(2) = toler(1)
  toler(3) = toler(1)*0.01d0
  toler(4) = toler(3)
  toler(5) = toler(3) 

  call get_input_rea(finp,'NUMERICAL','mitve',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  maxit(1) = rvoid(1)
  maxit(2) = maxit(1)

  call get_input_rea(finp,'NUMERICAL','mittu',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
  if(istat.lt.0) call runend(message)
  maxit(3) = rvoid(1)
  maxit(4) = maxit(3)

  
  call get_input_rea(finp,'NUMERICAL','damping',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
 ! if(istat.lt.0) call runend(message)
 if (istat.eq.0)  damping= rvoid(1)

  call get_input_rea(finp,'NUMERICAL','zdamping',rvoid,1,istat,message)
  if(istat.gt.0) call wriwar(message)
!  if(istat.lt.0) call runend(message)
  if(istat.eq.0)  z_damping = rvoid(1)
  !
  !*** Reads the THERMAL coupling block
  !
  ! Thermal model
  call get_input_int(finp,'THERMAL','thmod',ivoid,1,istat,message)
  if(istat.gt.0) then
     call wriwar(message)
  else if(istat.lt.0) then
     call wriwar(message)
  else if(ivoid(1).eq.1.or.ivoid(1).eq.4) then !index 4 is for retro-compatibility
     kfl_thmod = 1 
     if (ivoid(1).eq.4) print *,'readinp: thmod = 1,2,3 are deprecated, Sogachev model is 1'
     print *,'readinp: Thermal coupling activated, using Sogachev-Koblitz model.'     
  else if(ivoid(1).ne.0) then
     call runend('Thermal coupling model is incorrect or not implemented')
  end if

  ! Type of coupling
  if (kfl_thmod.ne.0) then
     kfl_thcou = .true. ! thermal coupling

     ! Time-varying surface temperature?
     call get_input_cha(finp,'THERMAL','trans',cvoid,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.&
          TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') &
          kfl_trtem =.true.
     
     if (kfl_trtem) then
        !Transient surface temperature cases
        call get_input_cha(finp,'THERMAL','case',cvoid,1,istat,message)
        if(istat.gt.0) then
           call wriwar(message)
        else if (TRIM(cvoid).eq.'gabls2'.or.TRIM(cvoid).eq.'GABLS2') then
             kfl_case = 1
        else if (TRIM(cvoid).eq.'cycle'.or.TRIM(cvoid).eq.'CYCLE') then
             kfl_case = 2
        else if (TRIM(cvoid).eq.'gabls3'.or.TRIM(cvoid).eq.'GABLS3') then
           kfl_case = 3
           print *,'readinp: Transient thermal case GABLS3 type'
        else if (TRIM(cvoid).eq.'gabls1'.or.TRIM(cvoid).eq.'GABLS1') then
           kfl_case = 4
        else
           kfl_case = 0
!           call runend('readinp: Transient THERMAL case is incorrect or not implemented')
        end if
        
        !GABLS3 case type
        if (kfl_case.eq.3) then
           call get_input_cha(finp,'THERMAL','tsfc',cvoid,1,istat,message)
           if(istat.gt.0) call wriwar(message)
           if(istat.lt.0) call runend(message)
           if (TRIM(cvoid).eq.'T2m_qw'.or.TRIM(cvoid).eq.'t2m_qw') &
                kfl_temp = 1
           if (TRIM(cvoid).eq.'tsk'.or.TRIM(cvoid).eq.'TSK') &
                kfl_temp = 2
           if (TRIM(cvoid).eq.'T2m'.or.TRIM(cvoid).eq.'t2m') &
                kfl_temp = 3

           call get_input_cha(finp,'THERMAL','th_adv',cvoid,1,istat,message)
           if(istat.gt.0) call wriwar(message)
           if(istat.lt.0) call runend(message)
           if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.&
                TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') &
                kfl_thadv =.true.

           call get_input_cha(finp,'THERMAL','mom_adv',cvoid,1,istat,message)
           if(istat.gt.0) call wriwar(message)
           if(istat.lt.0) call runend(message)
           if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.&
                TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') &
                kfl_momadv =.true.

           call get_input_cha(finp,'THERMAL','press_gr',cvoid,1,istat,message)
           if(istat.gt.0) call wriwar(message)
           if(istat.lt.0) call runend(message)
           if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.&
                TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') &
                kfl_pressgr =.true.

           call get_input_cha(finp,'THERMAL','nudging_mom',cvoid,1,istat,message)
           if(istat.gt.0) call wriwar(message)
           if(istat.eq.0) then
              if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.&
                   TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') &
                   kfl_mom_nudging =.true.
           end if

           call get_input_cha(finp,'THERMAL','nudging_temp',cvoid,1,istat,message)
           if(istat.gt.0) call wriwar(message)
           if(istat.eq.0) then
              if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.&
                   TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') &
                   kfl_tem_nudging =.true.
           end if

           call get_input_cha(finp,'THERMAL','meso_file',cvoid,1,istat,message)
           if(istat.gt.0) call wriwar(message)
           if(istat.lt.0) call runend(message)
           meso_file = TRIM(cvoid)
           print *, 'Meso file =', TRIM(meso_file)
        end if
     end if
     if (.not.kfl_trtem.or.kfl_case.eq.4.or.kfl_case.eq.0) then
        call get_input_rea(finp,'THERMAL','heatfl',rvoid,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        hflx0 = rvoid(1) ! heat flux over the ground
        if (abs(hflx0).gt.0.01) kfl_thcou = .true. 
        if (abs(hflx0).lt.0.01) kfl_thmod = 0 !stationary

        call get_input_rea(finp,'THERMAL','height',rvoid,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        ztsbl = rvoid(1) ! z top sbl
        if (kfl_thmod.ne.1) ztsbl =1.0d50

        call get_input_rea(finp,'THERMAL','molen',rvoid,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        lmoni = rvoid(1) !  Monin - Obukhov length (positive: stable, negative: unstable)

        call get_input_rea(finp,'THERMAL','tewall',rvoid,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        tewal = rvoid(1) 

        call get_input_rea(finp,'THERMAL','gradtop',rvoid,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        gradto = rvoid(1) 

        call get_input_rea(finp,'THERMAL','gradbot',rvoid,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        gradbo = rvoid(1) 

        call get_input_rea(finp,'THERMAL','zinfle',rvoid,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) call runend(message)
        ztemin = rvoid(1) 

        if (ztemin.lt.0.0.or.ztemin.gt.length) ztemin = length
        ! defines MO function parameters
        alpha_mo = 16.0d0  
        beta_mo  = 5.00d0
        expon_mo = - 0.25d0
     end if
  end if
  
  !
  !*** Reads the POSTPROCESS block
  !  
  if (kfl_trtem) then
     ! Postproces at some time steps
     call get_input_rea(finp,'POSTPROCESS','steps',rvoid,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     stepr = rvoid(1)
     iiter = 0

     !Tracking point (height)
     call get_input_rea(finp,'POSTPROCESS','cutheight',rvoid,mcuts,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)

     cutpr(1:mcuts) = rvoid(1:mcuts)
     do icuts =1, mcuts
        if (cutpr(icuts).le.dwall.or.cutpr(icuts).gt.length) &
             exit
     end do
     ncuts = icuts -1
     print *, 'ncuts', ncuts

     ! Output in netCDF
     call get_input_cha(finp,'POSTPROCESS','netCDF',cvoid,1,istat,message)
     if(istat.lt.0) then
        call wriwar(message)
     elseif (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
        kfl_netCDF =.true.
        print *,'readinp: netCDF output activated'
     end if

     if (kfl_netCDF) then
        !Output frequency in seconds (it should be a multiple of the time step)
        call get_input_int(finp,'POSTPROCESS','freq',ivoid,1,istat,message)
        if(istat.gt.0) call wriwar(message)
        if(istat.lt.0) then
           call wriwar(message)
           freq_netcdf = stepr  
           print *,'readinp: netCDF output frequency set to steps'
        end if
        if(istat.eq.0) freq_netcdf = ivoid(1)
!!$        if(mod(freq_netcdf,1.0d0/dtinv)/=0) then
!!$           call runend('readinp: netCDF output frequency is not a multiple of time step')
!!$        end if
        !Appends outputs to existing netCDF file
        if (restart_in) then
           call get_input_cha(finp,'POSTPROCESS','append',cvoid,1,istat,message)
           if(istat.gt.0) call wriwar(message)
           if(istat.lt.0) call runend(message)
           if (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes') then
              append_netCDF =.true.
              print *,'readinp: netCDF output append option activated'
           end if
        end if
     end if
  end if

  !
  !*** Reads the CANOPY block
  !
  call get_input_cha(finp,'CANOPY','canopy',cvoid,1,istat,message)
  if(istat.lt.0) then
     call wriwar(message)
     cvoid(1:2)='no'
  elseif (TRIM(cvoid).eq.'YES'.or.TRIM(cvoid).eq.'yes'.or.TRIM(cvoid).eq.'On'.or.TRIM(cvoid).eq.'on') then
     kfl_canop =.true.
  end if
     
  if (kfl_canop) then
     call get_input_rea(finp,'CANOPY','cd',rvoid,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     cdcan = rvoid(1) ! drag coeff

     call get_input_rea(finp,'CANOPY','LAD',rvoid,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     LAD = rvoid(1) ! Leaf area density

     call get_input_rea(finp,'CANOPY','height',rvoid,1,istat,message)
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     heica = rvoid(1) ! Height of the canopy 

     call get_input_cha(finp,'CANOPY','model',cvoid,1,istat,message)
     !     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) then
        call wriwar(message)
        call wriwar('Set Sogachev canopy model (default)')
     else if (TRIM(cvoid).eq.'SVENS') then
        kfl_canmo =1
        print *, 'CANOPY MODEL=SVENS'
     else if (TRIM(cvoid).eq.'SANZ') then
        kfl_canmo =2
        print *, 'CANOPY MODEL=SANZ'
     else if (TRIM(cvoid).eq.'GREEN') then
        kfl_canmo =3
        print *, 'CANOPY MODEL=GREEN'
     else if (TRIM(cvoid).eq.'LOPES') then
        kfl_canmo =4
        print *, 'CANOPY MODEL=LOPES'
     end if

     call get_input_cha(finp,'CANOPY','distribution',cvoid,1,istat,message)
     if(istat.lt.0) then
        call wriwar(message)
        call wriwar('Set uniform canopy distribution (default)')
     else if (TRIM(cvoid).eq.'LALIC') then ! from Lalic and mihailovic
        kfl_candi =1
     else if (TRIM(cvoid).eq.'TABLE') then  ! Table distr. (Lopes da Costa phd-thesis)
        kfl_candi =2      
     end if
     call get_input_rea(finp,'CANOPY','LAI',rvoid,1,istat,message)        
     if(istat.gt.0) call wriwar(message)
     if(istat.lt.0) call runend(message)
     LAI = rvoid(1) ! Leaf area index

     call get_input_cha(finp,'CANOPY','rad_file',cvoid,1,istat,message)
     if(istat.gt.0) call wriwar(message)
!     if(istat.lt.0) call runend(message)
     rad_file = TRIM(cvoid)
     print *, 'Radiation file =', TRIM(rad_file)

  end if

  return
end subroutine readinp
