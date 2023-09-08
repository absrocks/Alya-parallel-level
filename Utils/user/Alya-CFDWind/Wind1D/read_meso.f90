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
  implicit none

  integer(ip)                     :: ncID,VarID,DimID
  integer(ip)                     :: i
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
  character(len=s_name)           :: qw_name    = "wt"     !kinematic heat flux
  

  !
  ! Opens netCDF file
  !
  fmeso = TRIM(meso_file)
  if( nf90_open(TRIM(fmeso),NF90_NOWRITE, ncID) /= 0 ) &
       call runend('Error in nf90_open for file '//TRIM(fmeso))
  !
  ! Get times 
  !
  if( nf90_inq_dimid(ncID,time_dim,DimID) /= 0 ) &
       call runend('Error in nf90_inq_dimid '//TRIM(time_dim))
  if( nf90_inquire_dimension(ncID,DimID,time_dim,nt) /= 0 ) &
       call runend('Error in nf90_inquire_dimension '//TRIM(time_dim))
  !
  allocate(times(nt))
  allocate(seconds(nt))
  !
  if( nf90_inq_varid(ncID,time_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(time_name))
  if( nf90_get_var  (ncID,VarID,times,start=(/1/),count=(/nt/)) /= 0) &
       call runend('Error reading variable time')
  !
  ! Time manipulation
  ! comes like days since year 0 -> seconds since 30/06/2006 at 12:00  
  seconds(1) = 0_ip
  do i=2,nt
     seconds(i) = nint((times(i)*24.0_rp*3600.0_rp)-(times(1)*24.0_rp*3600.0_rp))
  end do
  !
  ! Get location
  !
  if( nf90_inq_varid(ncID,lon_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(lon_name))
  if( nf90_get_var  (ncID,VarID,lon) /= 0) &
       call runend('Error reading variable lon')
  !
    if( nf90_inq_varid(ncID,lat_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(lat_name))
  if( nf90_get_var  (ncID,VarID,lat) /= 0) &
       call runend('Error reading variable lat')
  !
  ! Get vertical levels
  !
  if( nf90_inq_dimid(ncID,z_dim,DimID) /= 0 ) &
       call runend('Error in nf90_inq_dimid '//TRIM(z_dim))
  if( nf90_inquire_dimension(ncID,DimID,z_dim,nz) /= 0 ) &
       call runend('Error in nf90_inquire_dimension '//TRIM(z_dim))
  !
  allocate(height(nz))
  !
  if( nf90_inq_varid(ncID,z_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(z_name))
  if( nf90_get_var  (ncID,VarID,height,start=(/1/),count=(/nz/)) /= 0) &
       call runend('Error reading variable z')
  !
  ! Get coriolis term
  !
  if( nf90_inq_varid(ncID,fc_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(fc_name))
  if( nf90_get_var  (ncID,VarID,fc) /= 0) &
       call runend('Error reading variable fc')
  fcori = fc
  write(*,*) "-Coriolis force taken from netCDF file-"
  !
  ! Get geostrophic wind
  !
  allocate(ugeo(nz,nt))
  allocate(vgeo(nz,nt))
  !
  if( nf90_inq_varid(ncID,ug_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(ug_name))
  if( nf90_get_var  (ncID,VarID,ugeo,start=(/1,1/),count=(/nz,nt/)) /= 0) &
       call runend('Error reading variable ugeo')
  !
  if( nf90_inq_varid(ncID,vg_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(vg_name))
  if( nf90_get_var  (ncID,VarID,vgeo,start=(/1,1/),count=(/nz,nt/)) /= 0) &
       call runend('Error reading variable vgeo')
  !
  !Get advections
  !
  allocate(uadv(nz,nt))
  allocate(vadv(nz,nt))
  allocate(thadv(nz,nt))
  !
  if( nf90_inq_varid(ncID,uadv_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(uadv_name))
  if( nf90_get_var  (ncID,VarID,uadv,start=(/1,1/),count=(/nz,nt/)) /= 0) &
       call runend('Error reading variable uadv')
 
  !
  if( nf90_inq_varid(ncID,vadv_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(vadv_name))
  if( nf90_get_var  (ncID,VarID,vadv,start=(/1,1/),count=(/nz,nt/)) /= 0) &
       call runend('Error reading variable uadv')
  
  !
  if( nf90_inq_varid(ncID,thadv_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(thadv_name))
  if( nf90_get_var  (ncID,VarID,thadv,start=(/1,1/),count=(/nz,nt/)) /= 0) &
       call runend('Error reading variable thadv')
  !
  ! Get near Surface Temperatures 
  !
  allocate(t2(nt))
  allocate(tsk(nt))
  !
  if( nf90_inq_varid(ncID,t2_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(t2_name))
  if( nf90_get_var  (ncID,VarID,t2,start=(/1/),count=(/nt/)) /= 0) &
       call runend('Error reading variable t2')
  !
  if( nf90_inq_varid(ncID,tsk_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(tsk_name))
  if( nf90_get_var  (ncID,VarID,tsk,start=(/1/),count=(/nt/)) /= 0) &
       call runend('Error reading variable tsk')
  !
  ! Get intial state
  !
  allocate(uini(nz))
  allocate(vini(nz))
  allocate(thini(nz))
  allocate(ustini(nt))
  allocate(hflux(nt))
  !U velocity component
  if( nf90_inq_varid(ncID,u_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(u_name))
  if( nf90_get_var  (ncID,VarID,uini,start=(/1,1/),count=(/nz,1/)) /= 0) &
       call runend('Error reading variable u')
  !V velocity component
  if( nf90_inq_varid(ncID,v_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(v_name))
  if( nf90_get_var  (ncID,VarID,vini,start=(/1,1/),count=(/nz,1/)) /= 0) &
       call runend('Error reading variable v')
  !Potential Temperature
  if( nf90_inq_varid(ncID,th_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(th_name))
  if( nf90_get_var  (ncID,VarID,thini,start=(/1,1/),count=(/nz,1/)) /= 0) &
       call runend('Error reading variable th')
  if (height(1).lt.1.0d-4) then
       write(*,*) "-Potential temperature at 0m height subtituted by tsk or t2m-"
       if (kfl_temp.eq.1.or.kfl_temp.eq.3) thini(1) = t2(2)
       if (kfl_temp.eq.2) thini(1) = tsk(2)
  end if
  !Friction Velocity
  if( nf90_inq_varid(ncID,ust_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(ust_name))
  if( nf90_get_var  (ncID,VarID,ustini) /= 0) &
       call runend('Error reading variable ust')
  !Kinematic heat flux
  if( nf90_inq_varid(ncID,qw_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(qw_name))
  if( nf90_get_var  (ncID,VarID,hflux) /= 0) &
       call runend('Error reading variable qw')
  hflux(1) = hflux(2) !hflux(1) was wrong in the netCDF file
  !
  ! Get velocity for top Boundary, Potential Tempearture
  !
  allocate(u(nz,nt))
  allocate(v(nz,nt))
  allocate(th(nz,nt))
  !U velocity component
  if( nf90_inq_varid(ncID,u_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(u_name))
  if( nf90_get_var  (ncID,VarID,u,start=(/1,1/),count=(/nz,nt/)) /= 0) &
       call runend('Error reading variable u')
  !V velocity component
  if( nf90_inq_varid(ncID,v_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(v_name))
  if( nf90_get_var  (ncID,VarID,v,start=(/1,1/),count=(/nz,nt/)) /= 0) &
       call runend('Error reading variable v')
  !Potential Temperature
  if( nf90_inq_varid(ncID,th_name,VarID) /= 0) &
       call runend('Error getting VarID for '//TRIM(th_name))
  if( nf90_get_var  (ncID,VarID,th,start=(/1,1/),count=(/nz,nt/)) /= 0) &
       call runend('Error reading variable th')

  !
  ! Calculates l_max
  !
  i = 1
  do while (height(i).lt.length)
     i = i + 1
  end do
  l_max = 0.00027*sqrt(u(i,1)*u(i,1)+v(i,1)*v(i,1))/abs(fcori) !JBR new lmax

  return 
end subroutine read_meso




!!$subroutine read_measu
!!$!****************************************************************
!!$!*
!!$!* Reads the GABLS3 case type input files
!!$!* This subroutine reads the dimensions and variables inside the netCDF file.
!!$!* In section "Variables names", these are specified with the names from the test
!!$!* netCDF file. For each new case, one should check the variables names inside the
!!$!* mesoscale file and change the names bellow if they do not coincide.
!!$!* 
!!$!* To check the variables names and atributes in the netCDF file:
!!$!* ncdump -h "file_name"
!!$!* To inspect the values of certain variable:
!!$!* ncdump -v "variable_name" "file_name"
!!$!*
!!$!* Jordi Barcons - PhD Student - 11/07/2016
!!$!*
!!$!****************************************************************
!!$
!!$  use KindType
!!$  use def_master
!!$  use InpOut  
!!$  use netCDF
!!$  implicit none
!!$
!!$  integer(ip)                     :: ncID,VarID,DimID
!!$  integer(ip)                     :: i
!!$  character(len=s_name)           :: fmeso 
!!$  !netCDF file varialbles names
!!$  character(len=256)              :: time_dim   = 'Time [days]'  !time in days (for GABLS3 case is "time")
!!$  character(len=256)              :: z_dim      = 'Height [m]' !height levels  (for GABLS3 case is "z")
!!$  character(len=s_name)           :: time_name  = 'time'   !time in days since 01/01/0000
!!$  character(len=s_name)           :: lon_name   = 'lon'
!!$  character(len=s_name)           :: lat_name   = 'lat'
!!$  character(len=s_name)           :: z_name     = 'z'      !height levels
!!$  character(len=s_name)           :: fc_name    = 'fc'     !Coriolis term
!!$  character(len=s_name)           :: ug_name    = 'Ug'     !Ug = Press. grad. RHS term (divided by fc) V-component
!!$  character(len=s_name)           :: vg_name    = 'Vg'     !Vg = - Press. grad. RHS term (divided by fc) U-component
!!$  character(len=s_name)           :: uadv_name  = 'Uadv'   !Advective momentum RHS term (divided by fc) U-component
!!$  character(len=s_name)           :: vadv_name  = 'Vadv'   !Advective momentum RHS term (divided by fc) V-component
!!$  character(len=s_name)           :: thadv_name = 'Thadv'  !Potential temperature advective RHS term
!!$  character(len=s_name)           :: t2_name    = 'T2'     !2m temperature
!!$  character(len=s_name)           :: tsk_name   = 'TSK'    !Skin temperature (Surface temperature)
!!$  character(len=s_name)           :: u_name     = "U"      !U velocity component
!!$  character(len=s_name)           :: v_name     = "V"      !V velocity component
!!$  character(len=s_name)           :: th_name    = "Th"     !Potential temperature
!!$  character(len=s_name)           :: ust_name   = "ust"    !Friction velocity
!!$  character(len=s_name)           :: qw_name    = "wt"     !kinematic heat flux
!!$  
!!$
!!$  !
!!$  ! Opens netCDF file
!!$  !
!!$  fmeso = TRIM(meso_file)
!!$  if( nf90_open(TRIM(fmeso),NF90_NOWRITE, ncID) /= 0 ) &
!!$       call runend('Error in nf90_open for file '//TRIM(fmeso))
!!$  !
!!$  ! Get times 
!!$  !
!!$  if( nf90_inq_dimid(ncID,time_dim,DimID) /= 0 ) &
!!$       call runend('Error in nf90_inq_dimid '//TRIM(time_dim))
!!$  if( nf90_inquire_dimension(ncID,DimID,time_dim,nt) /= 0 ) &
!!$       call runend('Error in nf90_inquire_dimension '//TRIM(time_dim))
!!$  !
!!$  allocate(times(nt))
!!$  allocate(seconds(nt))
!!$  !
!!$  if( nf90_inq_varid(ncID,time_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(time_name))
!!$  if( nf90_get_var  (ncID,VarID,times,start=(/1/),count=(/nt/)) /= 0) &
!!$       call runend('Error reading variable time')
!!$  !
!!$  ! Time manipulation
!!$  ! comes like days since year 0 -> seconds since 30/06/2006 at 12:00  
!!$  seconds(1) = 0_ip
!!$  do i=2,nt
!!$     seconds(i) = nint((times(i)*24.0_rp*3600.0_rp)-(times(1)*24.0_rp*3600.0_rp))
!!$  end do
!!$  !
!!$  ! Get location
!!$  !
!!$  if( nf90_inq_varid(ncID,lon_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(lon_name))
!!$  if( nf90_get_var  (ncID,VarID,lon) /= 0) &
!!$       call runend('Error reading variable lon')
!!$  !
!!$    if( nf90_inq_varid(ncID,lat_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(lat_name))
!!$  if( nf90_get_var  (ncID,VarID,lat) /= 0) &
!!$       call runend('Error reading variable lat')
!!$  !
!!$  ! Get vertical levels
!!$  !
!!$  if( nf90_inq_dimid(ncID,z_dim,DimID) /= 0 ) &
!!$       call runend('Error in nf90_inq_dimid '//TRIM(z_dim))
!!$  if( nf90_inquire_dimension(ncID,DimID,z_dim,nz) /= 0 ) &
!!$       call runend('Error in nf90_inquire_dimension '//TRIM(z_dim))
!!$  !
!!$  allocate(height(nz))
!!$  !
!!$  if( nf90_inq_varid(ncID,z_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(z_name))
!!$  if( nf90_get_var  (ncID,VarID,height,start=(/1/),count=(/nz/)) /= 0) &
!!$       call runend('Error reading variable z')
!!$  !
!!$  ! Get coriolis term
!!$  !
!!$  if( nf90_inq_varid(ncID,fc_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(fc_name))
!!$  if( nf90_get_var  (ncID,VarID,fc) /= 0) &
!!$       call runend('Error reading variable fc')
!!$  fcori = fc
!!$  write(*,*) "-Coriolis force taken from netCDF file-"
!!$  !
!!$  ! Get geostrophic wind
!!$  !
!!$  allocate(ugeo(nz,nt))
!!$  allocate(vgeo(nz,nt))
!!$  !
!!$  if( nf90_inq_varid(ncID,ug_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(ug_name))
!!$  if( nf90_get_var  (ncID,VarID,ugeo,start=(/1,1/),count=(/nz,nt/)) /= 0) &
!!$       call runend('Error reading variable ugeo')
!!$  !
!!$  if( nf90_inq_varid(ncID,vg_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(vg_name))
!!$  if( nf90_get_var  (ncID,VarID,vgeo,start=(/1,1/),count=(/nz,nt/)) /= 0) &
!!$       call runend('Error reading variable vgeo')
!!$  !
!!$  !Get advections
!!$  !
!!$  allocate(uadv(nz,nt))
!!$  allocate(vadv(nz,nt))
!!$  allocate(thadv(nz,nt))
!!$  !
!!$  if( nf90_inq_varid(ncID,uadv_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(uadv_name))
!!$  if( nf90_get_var  (ncID,VarID,uadv,start=(/1,1/),count=(/nz,nt/)) /= 0) &
!!$       call runend('Error reading variable uadv')
!!$ 
!!$  !
!!$  if( nf90_inq_varid(ncID,vadv_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(vadv_name))
!!$  if( nf90_get_var  (ncID,VarID,vadv,start=(/1,1/),count=(/nz,nt/)) /= 0) &
!!$       call runend('Error reading variable uadv')
!!$  
!!$  !
!!$  if( nf90_inq_varid(ncID,thadv_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(thadv_name))
!!$  if( nf90_get_var  (ncID,VarID,thadv,start=(/1,1/),count=(/nz,nt/)) /= 0) &
!!$       call runend('Error reading variable thadv')
!!$  !
!!$  ! Get near Surface Temperatures 
!!$  !
!!$  allocate(t2(nt))
!!$  allocate(tsk(nt))
!!$  !
!!$  if( nf90_inq_varid(ncID,t2_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(t2_name))
!!$  if( nf90_get_var  (ncID,VarID,t2,start=(/1/),count=(/nt/)) /= 0) &
!!$       call runend('Error reading variable t2')
!!$  !
!!$  if( nf90_inq_varid(ncID,tsk_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(tsk_name))
!!$  if( nf90_get_var  (ncID,VarID,tsk,start=(/1/),count=(/nt/)) /= 0) &
!!$       call runend('Error reading variable tsk')
!!$  !
!!$  ! Get intial state
!!$  !
!!$  allocate(uini(nz))
!!$  allocate(vini(nz))
!!$  allocate(thini(nz))
!!$  allocate(ustini(nt))
!!$  allocate(hflux(nt))
!!$  !U velocity component
!!$  if( nf90_inq_varid(ncID,u_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(u_name))
!!$  if( nf90_get_var  (ncID,VarID,uini,start=(/1,1/),count=(/nz,1/)) /= 0) &
!!$       call runend('Error reading variable u')
!!$  !V velocity component
!!$  if( nf90_inq_varid(ncID,v_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(v_name))
!!$  if( nf90_get_var  (ncID,VarID,vini,start=(/1,1/),count=(/nz,1/)) /= 0) &
!!$       call runend('Error reading variable v')
!!$  !Potential Temperature
!!$  if( nf90_inq_varid(ncID,th_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(th_name))
!!$  if( nf90_get_var  (ncID,VarID,thini,start=(/1,1/),count=(/nz,1/)) /= 0) &
!!$       call runend('Error reading variable th')
!!$  if (height(1).lt.1.0d-4) then
!!$       write(*,*) "-Potential temperature at 0m height subtituted by tsk or t2m-"
!!$       if (kfl_temp.eq.1.or.kfl_temp.eq.3) thini(1) = t2(2)
!!$       if (kfl_temp.eq.2) thini(1) = tsk(2)
!!$  end if
!!$  !Friction Velocity
!!$  if( nf90_inq_varid(ncID,ust_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(ust_name))
!!$  if( nf90_get_var  (ncID,VarID,ustini) /= 0) &
!!$       call runend('Error reading variable ust')
!!$  !Kinematic heat flux
!!$  if( nf90_inq_varid(ncID,qw_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(qw_name))
!!$  if( nf90_get_var  (ncID,VarID,hflux) /= 0) &
!!$       call runend('Error reading variable qw')
!!$  hflux(1) = hflux(2) !hflux(1) was wrong in the netCDF file
!!$  !
!!$  ! Get velocity for top Boundary, Potential Tempearture
!!$  !
!!$  allocate(u(nz,nt))
!!$  allocate(v(nz,nt))
!!$  allocate(th(nz,nt))
!!$  !U velocity component
!!$  if( nf90_inq_varid(ncID,u_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(u_name))
!!$  if( nf90_get_var  (ncID,VarID,u,start=(/1,1/),count=(/nz,nt/)) /= 0) &
!!$       call runend('Error reading variable u')
!!$  !V velocity component
!!$  if( nf90_inq_varid(ncID,v_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(v_name))
!!$  if( nf90_get_var  (ncID,VarID,v,start=(/1,1/),count=(/nz,nt/)) /= 0) &
!!$       call runend('Error reading variable v')
!!$  !Potential Temperature
!!$  if( nf90_inq_varid(ncID,th_name,VarID) /= 0) &
!!$       call runend('Error getting VarID for '//TRIM(th_name))
!!$  if( nf90_get_var  (ncID,VarID,th,start=(/1,1/),count=(/nz,nt/)) /= 0) &
!!$       call runend('Error reading variable th')
!!$
!!$  !
!!$  ! Calculates l_max
!!$  !
!!$  i = 1
!!$  do while (height(i).lt.length)
!!$     i = i + 1
!!$  end do
!!$  l_max = 0.00027*sqrt(u(i,1)*u(i,1)+v(i,1)*v(i,1))/abs(fcori) !JBR new lmax
!!$
!!$  return 
!!$end subroutine read_measu
