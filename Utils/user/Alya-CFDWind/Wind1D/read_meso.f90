Subroutine Read_meso
!****************************************************************
!*
!* Reads the GABLS3 case type input files
!* This subroutine reads the dimensions and variables inside the netCDF file.
!* In section "Variables names", these are specified with the names from the test
!* netCDF file. For each new case, one should check the variables names inside the
!* mesoscale file and change the names bellow if they do not coincide.
!* 
!* To check the variables names and atributes in the netCDF file:
!* ncdump -h "file_name"
!* To inspect the values of certain variable:
!* ncdump -v "variable_name" "file_name"
!*
!* Jordi Barcons - PhD Student - 11/07/2016
!*
!****************************************************************

  use KindType
  use def_master
  use InpOut  
  use netCDF
  use mod_usenetcdf
 
  implicit none

  integer(ip)                     :: ncID,VarID,DimID
  integer(ip)                     :: i,k 
  character(len=s_name)           :: fmeso 
  !netCDF file varialbles names
  character(len=256)              :: time_dim   = 'Time [days]'  !time in days (for GABLS3 case is "time")
  character(len=256)              :: z_dim      = 'Height [m]' !height levels  (for GABLS3 case is "z")
  character(len=s_name)           :: time_name  = 'time'   !time in days since 01/01/0000
  character(len=s_name)           :: lon_name   = 'lon'
  character(len=s_name)           :: lat_name   = 'lat'
  character(len=s_name)           :: z_name     = 'z'      !height levels
  character(len=s_name)           :: fc_name    = 'fc'     !Coriolis term
  character(len=s_name)           :: ug_name    = 'Ug'     !Ug = Press. grad. RHS term (divided by fc) V-component
  character(len=s_name)           :: vg_name    = 'Vg'     !Vg = - Press. grad. RHS term (divided by fc) U-component
  character(len=s_name)           :: uadv_name  = 'Uadv'   !Advective momentum RHS term (divided by fc) U-component
  character(len=s_name)           :: vadv_name  = 'Vadv'   !Advective momentum RHS term (divided by fc) V-component
  character(len=s_name)           :: thadv_name = 'Thadv'  !Potential temperature advective RHS term
  character(len=s_name)           :: t2_name    = 'T2'     !2m temperature
  character(len=s_name)           :: tsk_name   = 'TSK'    !Skin temperature (Surface temperature)
  character(len=s_name)           :: u_name     = "U"      !U velocity component
  character(len=s_name)           :: v_name     = "V"      !V velocity component
  character(len=s_name)           :: th_name    = "Th"     !Potential temperature
  character(len=s_name)           :: ust_name   = "ust"    !Friction velocity
  character(len=s_name)           :: wt_name    = "wt"     !kinematic heat flux
  character(len=s_name)           :: w_name, hf_name, L_name, time_units='Days'
  real (rp )                      :: facto
  real(rp)                        :: aux1, aux2, days, twall
  
  if (name=='GABLS3') then
     kfl_read_hf=.false.
     kfl_read_fc =.true.
     kfl_read_ug =.true.
     kfl_read_uadv =.true.
     kfl_read_thadv =.true.
     kfl_read_u =.true.
     kfl_read_w   =.false.
     kfl_read_th =.true.
     kfl_read_t2 =.true.
     kfl_read_tsk =.true.
     kfl_read_ust =.true.
     kfl_read_wt =.true.
     kfl_read_lat =.true.
     kfl_read_hf =.false.
     kfl_read_L  =.false.
     
  else if (name(2:5)=='laiz') then ! Alaiz
     time_dim   = 'time'
     z_dim      = 'level'
     
     time_name  = 'time'
     time_units = 'Seconds'
     z_name     = 'height'      !height levels

     kfl_read_lat=.false.
     lon_name   = 'lon'
     lat_name   = 'lat'
     
     kfl_read_fc =.true.
     fc_name    = 'fc'     !Coriolis term

     kfl_read_u =.true.   
     u_name     = "U"      !U velocity component
     v_name     = "V"      !V velocity component
     kfl_read_w =.true.
     w_name     = "W"      ! W velocity component
     kfl_read_th =.true.
     th_name    = "POT"     !Potential temperature

     kfl_read_ug =.true.
     ug_name    = 'Ug'     !Ug = Press. grad. RHS term (divided by fc) V-comp
     vg_name    = 'Vg'     !Vg = - Press. grad. RHS term (divided by fc) U-comp

     kfl_read_uadv =.true.
     uadv_name  = 'UADV'   !Advective momentum RHS term (divided by fc) U-comp
     vadv_name  = 'VADV'   !Advective momentum RHS term (divided by fc) V-comp

     kfl_read_thadv =.true.
     thadv_name = 'POT_ADV'  !Potential temperature advective RHS term

   
     kfl_read_t2 =.true.
     t2_name    = 'T2'     !2m temperature

     
     kfl_read_tsk =.true.
     tsk_name   = 'TSK'    !Skin temperature (Surface temperature)

     kfl_read_ust =.true.
     ust_name   = "UST"    !Friction velocity

     
     kfl_read_wt =.false.
     wt_name   = "wt"     !kinematic turb heat flux
   
     kfl_read_hf =.true.
     hf_name   = "HFX"    !turb heat flux  , hf=rhocp*wt

     kfl_read_L =.false.   ! Monin-Obukhov length
     L_name    ="L"       !
  end if
  
  !
  ! Opens netCDF file
  !
  fmeso = TRIM(meso_file)
  if( nf90_open(TRIM(fmeso),NF90_NOWRITE, ncID) /= 0 ) &
       call runend('Error in nf90_open for file '//TRIM(fmeso))
  !
  ! Get Dimensions ( number of times and vertical levels
  !
  call read_dim_netcdf(ncID, time_dim, nt) ! time dimension nt
  call read_dim_netcdf(ncID, z_dim, nz) ! number of height levels

  
  allocate(times(nt))
  allocate(seconds(nt))
  allocate(height(nz))
  
  !
  
  !
  !
  ! Get times 
  !
  call readv_netcdf(ncID, time_name, times, nt)
  
  ! Time manipulation
  ! comes like days
  if (time_units(1:5)=='Days ') facto = 24.0_rp*3600.0_rp
  if (time_units(1:5)=='Hours') facto = 3600_rp
  if (time_units(1:5)=='Secon') facto = 1_rp
  do i=1,nt
     seconds(i) = nint(times(i)*facto)-(times(1)*facto)
  end do
    
  
  ! Get vertical levels
  !
  call readv_netcdf(ncID, z_name, height, nz)

  !
  ! Get location
  !
  if (kfl_read_lat) then
     call readv_netcdf(ncID, lon_name, lon)
     call readv_netcdf(ncID, lat_name, lat)
  end if
  !
  ! Get coriolis term
  !
  if (kfl_read_fc) then
     call readv_netcdf(ncID, fc_name, fc)
     
     fcori = fc
     write(*,*) "-Coriolis force read from netCDF file-, fc=", fcori
  end if
  !
  !Get advections
  !
  if (kfl_read_uadv) then
     allocate(uadv(nz,nt))
     allocate(vadv(nz,nt))
     call readv_netcdf(ncID, uadv_name, uadv, nz, nt)
     call readv_netcdf(ncID, vadv_name, vadv, nz, nt)
  end if
  if (kfl_read_thadv) then
     allocate(thadv(nz,nt))
     call readv_netcdf(ncID, thadv_name, thadv, nz, nt)
  end if
  !
  ! Get geostrophic wind
  !
  if (kfl_read_ug) then
     allocate(ugeo(nz,nt))
     allocate(vgeo(nz,nt))
     call readv_netcdf(ncID, ug_name, ugeo, nz, nt)
     call readv_netcdf(ncID, vg_name, vgeo, nz, nt)
  end if
  !
  ! Get near Surface Temperatures 
  !
  !
  if (kfl_read_t2) then
     allocate(t2(nt))
     call readv_netcdf(ncID, t2_name, t2, nt)
  end if
  if (kfl_read_tsk) then
     allocate(tsk(nt))
     call readv_netcdf(ncID, tsk_name, tsk, nt)
  end if
  if (kfl_read_ust) then
     allocate(ust_wrf(nt))
     call readv_netcdf(ncID, ust_name, ust_wrf, nt)
  end if
  
  !
  !Kinematic heat flux
  !
  if (kfl_read_wt) then
     allocate(hflux(nt))
     call readv_netcdf(ncID, wt_name, hflux,nt)
  else if (kfl_read_hf) then
     allocate(hflux(nt))
     call readv_netcdf(ncID, hf_name, hflux,nt)
     hflux(1:nt) = hflux(1:nt)/rhocp 
  end if
  if (kfl_read_L) then
     allocate(MO_WRF(nt))
     call readv_netcdf(ncID, L_name, MO_WRF,nt)   
  end if

  
  !hflux(1) = hflux(2) !hflux(1) was wrong in the netCDF file
  !
  ! Get velocity for top Boundary, Potential Tempearture
  !
  if (kfl_read_u) then
     allocate(uini(nz))
     allocate(vini(nz))
     allocate(u(nz,nt))
     allocate(v(nz,nt))
     call readv_netcdf(ncID, u_name, u, nz, nt)
     call readv_netcdf(ncID, v_name, v, nz, nt)
     if (kfl_read_w)   then
        allocate(w(nz,nt))
        call readv_netcdf(ncID, w_name, w, nz, nt)
     end if
     !
  ! Get intial state
  !
     uini(1:nz) = u(1:nz, 1)
     vini(1:nz) = v(1:nz, 1)
  !
  ! Calculates l_max
  !
     i = 1
     do while (height(i).lt.length)
        i = i + 1
     end do
     l_max = 0.00027*sqrt(u(i,1)*u(i,1)+v(i,1)*v(i,1))/abs(fcori) !JBR new lmax
     
  end if
  if (kfl_read_th) then
     allocate(th(nz,nt))
     allocate(thini(nz))
     call readv_netcdf(ncID, th_name, th, nz, nt)
     thini(1:nz) = th(1:nz, 1)     
     if (height(1).lt.1.0d-4) then
        write(*,*) "-Potential temperature at 0m height subtituted by tsk or t2m-"
        if (kfl_temp.eq.1.or.kfl_temp.eq.3) thini(1) = t2(2) ! temp a 2m
        if (kfl_temp.eq.2) thini(1) = tsk(2)   ! surface temper
     end if
  end if
  !
  !  **** POSTPROCESS TENDENCIES
  !
   open(90,FILE='meso_th_splot.txt') !,STATUS='new')
   open(91,FILE='meso_ug_splot.txt') !,STATUS='new')
   open(92,FILE='meso_ugadv_splot.txt') !,STATUS='new')
   open(93,FILE='meso_adv_splot.txt') !,STATUS='new')
   open(94,FILE='tendencies_splot.txt') !,STATUS='new')
   open(95,FILE='meso_vel_splot.txt') !,STATUS='new')
   open(96,FILE='meso_WRFtime_splot.txt') !,STATUS='new')
   write(90,'(a)') '# 1: time [days] 2: height 3: th(heights) '
   write(91,'(a)') '# 1: time [days] 2: height 3: ug(heights) '
   write(92,'(a)') '# 1: time [days] 2: height 3: ug(heights) '
   write(93,'(a)') '# 1: time [days] 2: height 3: ug(heights) '
   write(94,'(a)') '# 1: time [days] 2: height 3: advx  4: advy  5: gpx  6: gpy  '
   write(95,'(a)') '# 1: time [days] 2: height 3: velx  4: vely  5: velz   '
   write(96,'(a)') '# 1: time [days] 2: T1     3: T2    4: TSK   5: T2m   6: hflux '
   do i =1, nt
      do k =1, nz
         days =  seconds(i)/(3600.0_rp*24.0_rp)
         write(90,'(500(e15.7,2x))') days, height(k) ,th(k,i)
         write(91,'(500(e15.7,2x))') days, height(k) ,sqrt(ugeo(k,i)*ugeo(k,i)+vgeo(k,i)*vgeo(k,i))
         aux1 = -vgeo(k,i) + uadv(k,i)
         aux2 = ugeo(k,i) + vadv(k,i)
         
         write(92,'(500(e15.7,2x))') days, height(k) ,sqrt(aux1*aux1 + aux2*aux2)
  
         write(93,'(500(e15.7,2x))') days, height(k) ,sqrt(uadv(k,i)*uadv(k,i)+vadv(k,i)*vadv(k,i))
         
         write(94,'(500(e15.7,2x))') days, height(k),  uadv(k,i), vadv(k,i), -vgeo(k,i),  ugeo(k,i)
         write(95,'(500(e15.7,2x))') days, height(k),  u(k,i), v(k,i)
         
      end do
      write (90,*)
      write (91,*)
      write (92,*)
      write (93,*)
      write (94,*)
      write (95,*)
      Twall = T2(i) + hflux(i)/(0.41*ust_wrf(i))*log(1.0 + 2.0/0.03)
      
      write (96,'(500(e15.7,2x))') days, th(1,i), th(2,i), tsk(i), T2(i), hflux(i), Twall !, 1.0/MO_WRF(i)
   end do
   close (90)
   close (91)
   close (92)
   close (93)
   close (94)
   close (95)
   close (96)
  return 
end subroutine read_meso



Subroutine Read_horna_tower_2
!****************************************************************
!*
!* Reads the GABLS3 case type input files
!* This subroutine reads the dimensions and variables inside the netCDF file.
!* In section "Variables names", these are specified with the names from the test
!* netCDF file. For each new case, one should check the variables names inside the
!* mesoscale file and change the names bellow if they do not coincide.
!* 
!* To check the variables names and atributes in the netCDF file:
!* ncdump -h "file_name"
!* To inspect the values of certain variable:
!* ncdump -v "variable_name" "file_name"
!*
!* Jordi Barcons - PostDoc Researcher - 24/05/2019
!*
!****************************************************************

  use KindType
  use def_master
  use InpOut  
  use netCDF
  use mod_usenetcdf
  implicit none

  integer(ip)                     :: ncID,VarID,DimID
  integer(ip)                     :: i, k
  character(len=s_name)           :: fmeso , vname
  !netCDF file varialbles names
  character(len=256)              :: time_dim   = 'Time [days]'  !time in days (for GABLS3 case is "time")
  character(len=256)              :: z_dim      = 'Height [m]' !height levels  (for GABLS3 case is "z")
  character(len=s_name)           :: lon_name   = 'lon'
  character(len=s_name)           :: lat_name   = 'lat'
  character(len=s_name)           :: fc_name    = 'fc'   !Coriolis term
  character(len=s_name)           :: time_name  = 'time' !time in days since 01/01/0000
  character(len=s_name)           :: z_name     = 'z'    !height levels
  character(len=s_name)           :: U_name    = 'U'     !U velocity component
  character(len=s_name)           :: V_name    = 'V'     !V velocity component
  character(len=s_name)           :: TKE_name    = 'TKE'
  character(len=s_name)           :: SWDOWN_name = 'SWDOWN'  !net incoming shortwave radiation
  character(len=s_name)           :: GLW_name    = 'GLW'   !net longwave radiation
  character(len=s_name)           :: ust_name   = "ust"    !Friction velocity
  character(len=s_name)           :: t2_name    = 'T2'     !2m temperature
  character(len=s_name)           :: tsk_name   = 'TSK'    !Skin temperature (Surface temperature)
  character(len=s_name)           :: hfx_name   = 'HFX'    !Heat Flux
  character(len=s_name)           :: wt_name    = "wt"     !kinematic heat flux
  character(len=s_name)           :: th_name    = "Th"     !Potential temperature  
  real(rp), allocatable           :: times2(:), height2(:)
  real(rp)                        :: aux1, aux2
  !
  ! Opens netCDF file
  !

  meso_file = 'Hornamossen_files/west_case/Hornamossen_interpolated_to_tower.nc'
  fmeso = TRIM(meso_file)
  if( nf90_open(TRIM(fmeso),NF90_NOWRITE, ncID) /= 0 ) &
       call runend('Error in nf90_open for file '//TRIM(fmeso))

  !
  ! Get dimensions
  !
  call read_dim_netcdf(ncID, time_dim, nt) ! time dimension nt
  call read_dim_netcdf(ncID, z_dim, nz)    ! vertical dimension nz

  !
  allocate(times(nt))
  allocate(seconds(nt))
  allocate(height(nz))
  allocate(t2(nt))
  allocate(tsk(nt))

  allocate(uini(nz))
  allocate(vini(nz))
  allocate(thini(nz))
  allocate(ust_wrf(nt))
  allocate(hflux(nt))
  allocate(ust(nt))

  allocate(u(nz,nt))
  allocate(v(nz,nt))
  allocate(th(nz,nt))

  allocate(glw(nt))   
  allocate(swdown(nt))
  allocate(hfx(nt))
 
  !
  ! Get times 
  !
  call readv_netcdf(ncID, time_name, times, nt)
  
  do i=1,nt
     seconds(i) = nint((times(i)*24.0_rp*3600.0_rp))
  end do

  !
  ! Get vertical levels
  !
  call readv_netcdf(ncID, z_name, height, nz)
  
  !
  ! Get coriolis term
  !
  call readv_netcdf(ncID, fc_name, fc)
  !
  ! Get location
  !
  call readv_netcdf(ncID, lon_name, lon)
  call readv_netcdf(ncID, lat_name, lat)
 
  fcori = fc
  write(*,*) "-Coriolis force taken from netCDF file-"

  !
  !
  ! Get near Surface Temperatures 
  !
  call readv_netcdf(ncID, t2_name, t2, nt)
  call readv_netcdf(ncID, tsk_name, tsk, nt)
 
  !Friction Velocity
  !
  call readv_netcdf(ncID, ust_name, ust, nt)
  call readv_netcdf(ncID, wt_name, hflux, nt)

  !
  ! Get velocities Potential Tempearture
  !
  
  call readv_netcdf(ncID, u_name, u, nz, nt)  !U velocity component
  call readv_netcdf(ncID, v_name, v, nz, nt)  !V velocity component
  call readv_netcdf(ncID, th_name, th, nz, nt)  !Potential Temperature
  
  ! Get radiative heat flux above canopy
  !
  
  call readv_netcdf(ncID, glw_name, glw, nt)  ! longwave radiation 
  call readv_netcdf(ncID, swdown_name, swdown, nt) ! shortwave radiation
  call readv_netcdf(ncID, hfx_name, hfx, nt) ! hfx
  
 
  ! write
   open(90,FILE='meso_input.txt') !,STATUS='new')
   write(90,'(a)') '# 1: time [days] 2: time[s] 3: glw 4: swdowm 5 :hflx  6: wt 7: ust 8: T2  9: TSK'
   do i =1, nt
      write(90,'((e15.7,2x),x,i,x,20(e15.7,2x))')  times(i) , seconds(i), glw(i), swdown(i), &
           hfx(i), hflux(i), ust(i), T2(i), TSK(i), th(1, i)
   end do

   close(90)


   open(90,FILE='meso_height.txt') !,STATUS='new')
   write(90,'(a)') '# 1: height [m] 2: ux 3: uy 4 : Temp00  5: Temp14'
   
   do i =1, nz
      write(90,'(20(e15.7,2x))')  height(i), u(i,200), v(i,200), &
           th(i, 1), th(i, 85)
   end do
   
   close(90)

   

   if (nf90_close(ncID) ==0) continue 


   fmeso = 'Hornamossen_files/west_case/Hornamossen_area_average.nc'
   if( nf90_open(TRIM(fmeso),NF90_NOWRITE, ncID) /= 0 ) &
        call runend('Error in nf90_open for file '//TRIM(fmeso))
   

  !
  ! Get dimensions
  !
!  call read_dim_netcdf(ncID, time_dim, nt) ! time dimension nt
!  call read_dim_netcdf(ncID, z_dim, nz)    ! vertical dimension nz


  !
  ! Get times 
  !
   allocate(times2(nt))
   allocate(height2(nz))
 
   call readv_netcdf(ncID, time_name, times2, nt)
   call readv_netcdf(ncID, z_name, height2, nz)
  
   print *, 'va a comparar'
   do i =1, nt
      if(abs(times(i)-times2(i)).gt.1.0e-7) call runend('uhhhloco mal los tiempos')
   end do
   do i =1, nz
      if(abs(height(i)-height2(i)).gt.1.0e-7) call runend('uhhhloco mal lah alturahs')
   end do
   
   print *, 'comparo bien'

   allocate( uadv(nz, nt))      
   allocate( vadv(nz, nt))
   allocate( ugeo(nz, nt))
   allocate( vgeo(nz, nt))
   
   vname='Uadv'
   call readv_netcdf(ncID, vname, uadv, nz, nt)
   vname='Vadv'
   call readv_netcdf(ncID, vname, vadv, nz, nt)
   vname='Ug'
   call readv_netcdf(ncID, vname, ugeo, nz, nt)
   vname='Vg'
   call readv_netcdf(ncID, vname, vgeo, nz, nt)

   open(90,FILE='meso_th_all.txt') !,STATUS='new')
   write(90,'(a)') '# 1: time [days] 2: th(heights) '
   do i =1, nt
      write(90,'(500(e15.7,2x))') times(i), (th(k,i) ,k=1,nz)
   end do
   print *, 'termino bien'
   close(90)

   open(90,FILE='meso_th_splot.txt') !,STATUS='new')
   open(91,FILE='meso_ug_splot.txt') !,STATUS='new')
   open(92,FILE='meso_ugadv_splot.txt') !,STATUS='new')
   open(93,FILE='meso_adv_splot.txt') !,STATUS='new')
   open(94,FILE='tendencies_splot.txt') !,STATUS='new')
   write(90,'(a)') '# 1: time [days] 2: height 3: th(heights) '
   write(91,'(a)') '# 1: time [days] 2: height 3: ug(heights) '
   write(94,'(a)') '# 1: time [days] 2: height 3: advx  4: advy  5: gpx  6: gpy  '
  
   do i =1, nt
      do k =1, nz
         write(90,'(500(e15.7,2x))') times(i), height(k) ,th(k,i)
         write(91,'(500(e15.7,2x))') times(i), height(k) ,sqrt(ugeo(k,i)*ugeo(k,i)+vgeo(k,i)*vgeo(k,i))
         aux1 = -vgeo(k,i) + uadv(k,i)
         aux2 = ugeo(k,i) + vadv(k,i)
         
         write(92,'(500(e15.7,2x))') times(i), height(k) ,sqrt(aux1*aux1 + aux2*aux2)
  
         write(93,'(500(e15.7,2x))') times(i), height(k) ,sqrt(uadv(k,i)*uadv(k,i)+vadv(k,i)*vadv(k,i))
         
         write(94,'(500(e15.7,2x))') times(i), height(k),  uadv(k,i), vadv(k,i), -vgeo(k,i),  ugeo(k,i)
         
         
      end do
      write (90,*)
      write (91,*)
      write (92,*)
      write (93,*)
      write (94,*)
   end do
  
   
   print *, 'termino bien'
   close(90)
   close(91)
   close(92)
   close(93)
   close(94)
   open(90,FILE='meso_prgr40.txt') !,STATUS='new')
   do i =1, nt
      write(90,'(500(e15.7,2x))') times(i) ,sqrt(ugeo(5,i)*ugeo(5,i)+vgeo(5,i)*vgeo(5,i))
   end do
   close(90)
   if (nf90_close(ncID) ==0) continue

   return
end subroutine Read_horna_tower_2

