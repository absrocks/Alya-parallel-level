subroutine netCDF_define
  !****************************************************************
  !*
  !* Defines netCDF dimensions and variables
  !*
  !* Jordi Barcons - PhD Student - 14/07/2016
  !*
  !****************************************************************

  use KindType
  use def_master
  use def_netCDF
  use InpOut  
  use netCDF
  implicit none
  logical :: define = .true.

  !*************************************************!
  !****************** DEFINE MODE ******************!
  !*************************************************!
  fnc_out = "Wind1D_output.nc"
  inquire(FILE=trim(fnc_out),EXIST=there)

  if (there.and..not.restart_in) then
     write(rm,'("rm ",A)')trim(fnc_out)
     call system(rm)
     define = .true.
  else if (.not.there) then
     define = .true.
  else if (there.and.restart_in.and.append_netCDF) then
     define = .false.
  else if (there.and.restart_in.and..not.append_netCDF) then
     fnc_out = "Wind1D_output_2.nc"
     define = .true. 
  end if

  write (*,*) 'netCDF define=',define
  if (define) then
     !
     !*** Open netCDF file (define mode) and get ncID
     !
     if( nf90_create(TRIM(fnc_out),cmode=or(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid=ncID) /= 0 ) & 
          call runend('netCDF_define : Error in nf90_create') 
     !
     !**  Define dimensions 
     !
     if( nf90_def_dim(ncID, z_name, npoin, nz_ID ) /= 0 ) &
          call runend('netCDF_define : error in nf90_def_dim for z')
     if( nf90_def_dim(ncID, time_name,nf90_unlimited, nt_ID ) /= 0 ) & 
          call runend('netCDF_define : error in nf90_def_dim for time')
     if( nf90_def_dim(ncID, site_name, 1, site_ID ) /= 0 ) & 
          call runend('netCDF_define : error in nf90_def_dim for site')
     !
     !*** Define coordinates variables
     !                 
     ! Heights
     if( nf90_def_var(ncID, z_name ,NF90_FLOAT, (/nz_ID/), z_ID) /= 0 ) &
          call runend('netCDF_define : error in nf90_def_variable for variable '//TRIM(z_name))
     attr_desc  = 'Height'     
     attr_units = 'm' 
     attr_coor = 'Height'  
     if( nf90_put_att(ncID, z_ID, 'units', attr_units) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     if( nf90_put_att(ncID, z_ID, 'description', attr_desc) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att') 
     if( nf90_put_att(ncID, z_ID, '_CoordinateAxisType', attr_coor) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     !
     ! Time in seconds
     if( nf90_def_var(ncID, time_name, NF90_INT, (/nt_ID/), time_ID) /= 0 ) &
          call runend('netCDF_define : error in nf90_def_var for variable '//TRIM(time_name))  
     attr_desc  = 'Time in seconds'
     attr_units = 's' 
     attr_coor = 'Time' 
     if( nf90_put_att(ncID, time_ID, 'units', attr_units) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     if( nf90_put_att(ncID, time_ID, 'description', attr_desc) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     if( nf90_put_att(ncID, time_ID, '_CoordinateAxisType', attr_coor) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     !
     !*** Define column variables
     !
     ! U
     if( nf90_def_var(ncID, u_name, NF90_FLOAT, (/nz_ID,nt_ID/), u_ID) /= 0 ) &
          call runend('outnc_meso : error in nf90_def_var for variable '//TRIM(u_name))  
     attr_desc  = 'U velocity'
     attr_units = 'm s-1'  
     if( nf90_put_att(ncID, u_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, u_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! V
     if( nf90_def_var(ncID, v_name, NF90_FLOAT, (/nz_ID,nt_ID/), v_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(v_name))  
     attr_desc  = 'V velocity'
     attr_units = 'm s-1'  
     if( nf90_put_att(ncID, v_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, v_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! VEL
     if( nf90_def_var(ncID, vel_name, NF90_FLOAT, (/nz_ID,nt_ID/), vel_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(vel_name))  
     attr_desc  = 'velocity module'
     attr_units = 'm s-1'  
     if( nf90_put_att(ncID, vel_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, vel_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! DIR
     if( nf90_def_var(ncID, dir_name, NF90_FLOAT, (/nz_ID,nt_ID/), dir_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(dir_name))  
     attr_desc  = 'velocity meteorological direction'
     attr_units = 'degrees'  
     if( nf90_put_att(ncID, dir_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, dir_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! T
     if( nf90_def_var(ncID, th_name, NF90_FLOAT, (/nz_ID,nt_ID/), th_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(th_name))
     attr_desc  = 'Potential Temperature'
     attr_units = 'K'
     if( nf90_put_att(ncID, th_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, th_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! U advection
     if( nf90_def_var(ncID, uadv_name, NF90_FLOAT, (/nz_ID,nt_ID/), uadv_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(uadv_name))  
     attr_desc  = 'Advective momentum RHS term (divided by fc) U-component'
     attr_units = 'm s-1'  
     if( nf90_put_att(ncID, uadv_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, uadv_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! V advection
     if( nf90_def_var(ncID, vadv_name, NF90_FLOAT, (/nz_ID,nt_ID/), vadv_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(vadv_name))  
     attr_desc  = 'Advective momentum RHS term (divided by fc) V-component'
     attr_units = 'm s-1'  
     if( nf90_put_att(ncID, vadv_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, vadv_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! T advection
     if( nf90_def_var(ncID, thadv_name, NF90_FLOAT, (/nz_ID,nt_ID/), thadv_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(thadv_name))
     attr_desc  = 'Potential temperature advective RHS term'
     attr_units = 'K s-1'
     if( nf90_put_att(ncID, thadv_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, thadv_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! TKE
     if( nf90_def_var(ncID, tke_name, NF90_FLOAT, (/nz_ID,nt_ID/), tke_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(tke_name))  
     attr_desc  = 'Turbulent Kinetic Energy'
     attr_units = 'm2 s-2'  
     if( nf90_put_att(ncID, tke_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, tke_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! Dissipation
     if( nf90_def_var(ncID, eps_name, NF90_FLOAT, (/nz_ID,nt_ID/), eps_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(eps_name))  
     attr_desc  = 'Dissipation'
     attr_units = 'm2 s-3'  
     if( nf90_put_att(ncID, eps_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, eps_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! Turbulent Viscosity
     if( nf90_def_var(ncID, mut_name, NF90_FLOAT, (/nz_ID,nt_ID/), mut_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(mut_name))  
     attr_desc  = 'Turbulent Viscosity'
     attr_units = 'm2 s-1'  
     if( nf90_put_att(ncID, mut_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, mut_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! Mixing Lenght
     if( nf90_def_var(ncID, lm_name, NF90_FLOAT, (/nz_ID,nt_ID/), lm_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(lm_name))  
     attr_desc  = 'Mixing Lenght'
     attr_units = 'm'  
     if( nf90_put_att(ncID, lm_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, lm_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')  
     !
     ! Flux Richardson Number
     if( nf90_def_var(ncID, Rf_name, NF90_FLOAT, (/nz_ID,nt_ID/), Rf_ID) /= 0 ) &
          call runend('outnc_meso : error in nf90_def_var for variable '//TRIM(Rf_name))  
     attr_desc  = 'Flux Richardson Number'
     attr_units = ''  
     if( nf90_put_att(ncID, Rf_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, Rf_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! TKE Mechanical Production
     if( nf90_def_var(ncID, prdm_name, NF90_FLOAT, (/nz_ID,nt_ID/), prdm_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(prdm_name))  
     attr_desc  = 'TKE  Production'
     attr_units = 'kg m-1 s-3'  
     if( nf90_put_att(ncID, prdm_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, prdm_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! TKE Thermal Production
     if( nf90_def_var(ncID, prdt_name, NF90_FLOAT, (/nz_ID,nt_ID/), prdt_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(prdt_name))  
     attr_desc  = 'TKE Thermal Production'
     attr_units = 'kg m-1 s-3'  
     if( nf90_put_att(ncID, prdt_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, prdt_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     !
     ! Heat Flux
     if( nf90_def_var(ncID, hfx_name, NF90_FLOAT, (/nz_ID,nt_ID/), hfx_ID) /= 0 ) &
          call runend('outnc_meso : error in nf90_def_var for variable '//TRIM(hfx_name))  
     attr_desc  = 'Heat Flux'
     attr_units = 'kg s-3'  
     if( nf90_put_att(ncID, hfx_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_meso : error in nf90_put_att')
     if( nf90_put_att(ncID, hfx_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_meso : error in nf90_put_att')
     !
     ! Shear stress
     if( nf90_def_var(ncID, strs_name, NF90_FLOAT, (/nz_ID,nt_ID/),strs_ID) /= 0 ) &
          call runend('outnc_define : error in nf90_def_var for variable '//TRIM(strs_name))  
     attr_desc  = 'Shear stress'
     attr_units = 'kg m-1 s-2'  
     if( nf90_put_att(ncID, strs_ID, 'units', attr_units) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')
     if( nf90_put_att(ncID, strs_ID, 'description', attr_desc) /= 0 ) &
          call runend('outnc_define : error in nf90_put_att')

     !
     !*** Define other variables
     !
     ! Friction Velocity
     if( nf90_def_var(ncID, ust_name, NF90_FLOAT, (/nt_ID/), ust_ID) /= 0 ) &
          call runend('netCDF_define : error in nf90_def_var for variable '//TRIM(ust_name))  
     attr_desc  = 'Surface friction velocity'
     attr_units = 'm s-1' 
     if( nf90_put_att(ncID, ust_ID, 'units', attr_units) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     if( nf90_put_att(ncID, ust_ID, 'description', attr_desc) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     !
     ! Kinematic Upward sensible heat flux
     if( nf90_def_var(ncID, sflux_name, NF90_FLOAT, (/nt_ID/), sflux_ID) /= 0 ) &
          call runend('netCDF_define : error in nf90_def_var for variable '//TRIM(sflux_name))  
     attr_desc  = 'Surface kinematic sensible heat flux wt = HFX/(rho*Cp)'
     attr_units = 'K m s-1' 
     if( nf90_put_att(ncID, sflux_ID, 'units', attr_units) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     if( nf90_put_att(ncID, sflux_ID, 'description', attr_desc) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     !
     ! Obukhov Length
     if( nf90_def_var(ncID, LMO_name, NF90_FLOAT, (/nt_ID/), LMO_ID) /= 0 ) &
          call runend('netCDF_define : error in nf90_def_var for variable '//TRIM(LMO_name))  
     attr_desc  = 'Obukhov Length'
     attr_units = 'm' 
     if( nf90_put_att(ncID, LMO_ID, 'units', attr_units) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     if( nf90_put_att(ncID, LMO_ID, 'description', attr_desc) /= 0 ) &
          call runend('netCDF_define : error in nf90_put_att')
     !
     ! 2m Temperature
     ! if( nf90_def_var(ncID, T2_name, NF90_FLOAT, (/nt_ID/), T2_ID) /= 0 ) &
     !      call runend('netCDF_define : error in nf90_def_var for variable '//TRIM(T2_name))  
     ! attr_desc  = '2-m temperature'
     ! attr_units = 'K' 
     ! if( nf90_put_att(ncID, T2_ID, 'units', attr_units) /= 0 ) &
     !      call runend('netCDF_define : error in nf90_put_att')
     ! if( nf90_put_att(ncID, T2_ID, 'description', attr_desc) /= 0 ) &
     !      call runend('netCDF_define : error in nf90_put_att')
     !
     !*** Leave the define mode
     !
     if( nf90_enddef(ncID) /= 0 ) & 
          call runend('netCDF_define: error in nf90_enddef')

     !FINISHES DEFINE
  else
     !
     !*** Open netCDF file (write mode)
     !
     if( nf90_open(TRIM(fnc_out),NF90_WRITE, ncID) /= 0 ) & 
          call runend('netCDF_open : Error in nf90_open')
     !
     ! Inquire time dimension of the file
     if( nf90_inq_dimid(ncID, time_name, nt_ID) /= 0 ) & 
          call runend('netCDF_define : Error in time dimension')
     if( nf90_inquire_dimension(ncID, nt_ID, time_name, nt_netCDF) /= 0 ) & 
          call runend('netCDF_define : Error in time dimension')
     write (*,*) 'netCDF time dimension =',nt_netCDF
     

  end if

  !
  ! Closes the file
  if( nf90_close(ncID) /= 0) call runend('netCDF_define : Error in nf90_close')
  
  call netCDF_allocate
  
  return
end subroutine netCDF_define
  

subroutine netCDF_output
  !****************************************************************
  !*
  !* Writes netCDF output file
  !*
  !* Jordi Barcons - PhD Student - 14/07/2016
  !*
  !****************************************************************

  use KindType
  use def_master
  use def_netCDF
  use InpOut  
  use netCDF
  implicit none

  
  !*********************************************************!
  !****************** FILL VARIABLES MODE ******************!
  !*********************************************************!
  !
  nt_netCDF = nt_netCDF + 1
  write(*,*) 'netCDF_output: nt_netCDF,ctime=',nt_netCDF,ctime

  
  !
  !*** Open netCDF file (write mode)
  !
  if( nf90_open(TRIM(fnc_out),NF90_WRITE, ncID) /= 0 ) & 
       call runend('netCDF_open : Error in nf90_open')

  !*** Coordinates
  !
  if(.not.append_netCDF) then
  if( nf90_put_var(ncID, z_ID, coord(1:npoin)) /= 0 )  &
       call runend('setnc: error in nf90_put_var for varialbe z')
  end if
  !
  !
  !*** Writing variables
  ! Time
  istat = nf90_inq_varid(ncID,time_name,time_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(time_name))
  istat = nf90_put_var(ncID,time_ID, time_out,start=(/nt_netCDF/))
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(time_name))  
  ! U
  istat = nf90_inq_varid(ncID,u_name,u_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(u_name))
  istat = nf90_put_var(ncID,u_ID,  u_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(u_name))
  !
  ! V
  istat = nf90_inq_varid(ncID,v_name,v_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(v_name))
  istat = nf90_put_var(ncID,v_ID, v_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(v_name))
  !
  ! VEL
  istat = nf90_inq_varid(ncID,vel_name,vel_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(vel_name))
  istat = nf90_put_var(ncID,vel_ID, vel_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(vel_name))
  !
  ! DIR
  istat = nf90_inq_varid(ncID,dir_name,dir_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(dir_name))
  istat = nf90_put_var(ncID,dir_ID, dir_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(v_name))  
  !
  ! Th
  istat = nf90_inq_varid(ncID,th_name,th_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(th_name))
  istat = nf90_put_var(ncID,th_ID, th_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(th_name))
  !
  ! Uadv  
  istat = nf90_inq_varid(ncID,uadv_name,uadv_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(uadv_name))
  istat = nf90_put_var(ncID,uadv_ID, uadv_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(uadv_name))
  !
  ! Vadv
  istat = nf90_inq_varid(ncID,vadv_name,vadv_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(vadv_name))
  istat = nf90_put_var(ncID,vadv_ID, vadv_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(vadv_name))
  !
  ! Thadv
  istat = nf90_inq_varid(ncID,thadv_name,thadv_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(thadv_name))
  istat = nf90_put_var(ncID,thadv_ID, thadv_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(thadv_name))
  !
  ! TKE
  istat = nf90_inq_varid(ncID,tke_name,tke_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(tke_name))
  istat = nf90_put_var(ncID,tke_ID, tke_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(tke_name))
  !
  ! EPS
  istat = nf90_inq_varid(ncID,eps_name,eps_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(eps_name))
  istat = nf90_put_var(ncID,eps_ID, eps_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(eps_name))
  !
  ! MUT
  istat = nf90_inq_varid(ncID,mut_name,mut_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(mut_name))
  istat = nf90_put_var(ncID,mut_ID, mut_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(mut_name))
  !
  ! LM
  istat = nf90_inq_varid(ncID,lm_name,lm_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(lm_name))
  istat = nf90_put_var(ncID,lm_ID, lm_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(lm_name))
  !
  ! Rf
  istat = nf90_inq_varid(ncID,Rf_name,Rf_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(Rf_name))
  istat = nf90_put_var(ncID,Rf_ID, Rf_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(Rf_name))
  !
  ! PRDM
  istat = nf90_inq_varid(ncID,prdm_name,prdm_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(prdm_name))
  istat = nf90_put_var(ncID,prdm_ID, prdm_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(prdm_name))
  !
  ! PRDT
  istat = nf90_inq_varid(ncID,prdt_name,prdt_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(prdt_name))
  istat = nf90_put_var(ncID,prdt_ID, prdt_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(prdt_name))
  !
  ! HFX
  istat = nf90_inq_varid(ncID,hfx_name,hfx_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(hfx_name))
  istat = nf90_put_var(ncID,hfx_ID, hfx_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(hfx_name))
  !
  !STRS
  istat = nf90_inq_varid(ncID,strs_name,strs_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(strs_name))
  istat = nf90_put_var(ncID,strs_ID, strs_out(1:npoin), start=(/1,nt_netCDF/), &
       count=(/npoin,1/) )
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(strs_name))
  !
  !*** Write other variables
  !
  ! UST
  istat = nf90_inq_varid(ncID,ust_name,ust_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(ust_name))
  istat = nf90_put_var(ncID,ust_ID, ust_out, start=(/nt_netCDF/))
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(ust_name))  
  !
  ! SFLUX
  istat = nf90_inq_varid(ncID,sflux_name,sflux_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(sflux_name))
  istat = nf90_put_var(ncID,sflux_ID, sflux_out, start=(/nt_netCDF/))
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(sflux_name))
  !
  ! LMO
  istat = nf90_inq_varid(ncID,LMO_name,LMO_ID)
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(LMO_name))
  istat = nf90_put_var(ncID,LMO_ID, LMO_out, start=(/nt_netCDF/))
  if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(LMO_name))
  !
  ! T2
  ! istat = nf90_inq_varid(ncID,T2_name,T2_ID)
  ! if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(T2_name))
  ! istat = nf90_put_var(ncID,T2_ID, T2_out, start=(/nt_netCDF/))
  ! if(istat /= 0) call wriwar(nf90_strerror(istat)//'writing '//TRIM(T2_name))

  !
  ! Closes the file
  if( nf90_close(ncID) /= 0) call runend('netCDF_output : Error in nf90_close') 

  return
end subroutine netCDF_output


subroutine netCDF_Postpr
  ! Filling netCDF output variables
  use def_master
  use def_netCDF
  implicit none
  real(8)             :: tstar,  gpnut, lenmo
  real(8)             :: vel, mut,mixle,Richf,prodm,prodt,htflx,stres
  integer(4 )         :: ipoin

  gpnut =  cmu*keyva(1,1)*keyva(1,1)/epsil(1,1)
  ! qwall positive in stable atm, when floor is cold (heat flux going out)
  hflx0 =  rhocp*gpnut*(tempe(2,1)- tempe(1,1))/(coord(2)-coord(1))/sigte
  tstar =  hflx0/(rhocp*ustar)
  if (abs(tstar).lt.1.0e-8) tstar = 1.0e-8
  lenmo = ustar*ustar*teref/(kar*gravi*tstar)
  !
  time_out        = ctime
  ust_out         = ustar
  sflux_out       = hflx0/rhocp
  LMO_out         = lenmo
  u_out  (1:npoin) = veloc(1:npoin,1,1)
  v_out  (1:npoin) = veloc(1:npoin,2,1)
  vel_out(1:npoin) = sqrt(veloc(1:npoin,1,1)**2+veloc(1:npoin,2,1)**2)
  !
  do ipoin=1,npoin
     dir_out(ipoin) = (180.0d0/pi)*(atan2(-veloc(ipoin,1,1),-veloc(ipoin,2,1)))
     if (dir_out(ipoin).lt.0.0d0) dir_out(ipoin) = dir_out(ipoin) + 360.0d0
     call Postpr_ipoin(vel,mut,mixle,Richf,prodm,prodt,htflx,stres, ipoin)
     mut_out(ipoin)  = mut
     lm_out(ipoin)   = mixle
     Rf_out(ipoin)   = Richf
     prdm_out(ipoin) = prodm
     prdt_out(ipoin) = prodt
     hfx_out(ipoin)  = htflx
     strs_out(ipoin) = stres
  end do
  !
  th_out   (1:npoin) = tempe(1:npoin,1)
  !
  if (kfl_case.eq.3) then
     uadv_out (1:npoin) = u_adv(1:npoin)
     vadv_out (1:npoin) = v_adv(1:npoin)
     thadv_out(1:npoin) = tempe_adv(1:npoin)
  else
     uadv_out (1:npoin) = 0.0d0
     vadv_out (1:npoin) = 0.0d0
     thadv_out(1:npoin) = 0.0d0
  end if
  !
  tke_out  (1:npoin) = keyva(1:npoin,1)
  eps_out  (1:npoin) = epsil(1:npoin,1)
 
  call netCDF_output

end subroutine netCDF_Postpr


subroutine netCDF_allocate
  !
  ! Allocates netCDF output variables
  !
  use def_master
  use def_netCDF
  implicit none

  allocate (time_out)
  allocate (u_out     (npoin))
  allocate (v_out     (npoin))
  allocate (vel_out   (npoin))
  allocate (dir_out   (npoin))
  allocate (th_out    (npoin))
  allocate (uadv_out  (npoin))
  allocate (vadv_out  (npoin))
  allocate (thadv_out (npoin))
  allocate (tke_out   (npoin))
  allocate (eps_out   (npoin))
  allocate (mut_out   (npoin))
  allocate (lm_out    (npoin))
  allocate (Rf_out    (npoin))
  allocate (prdm_out  (npoin))
  allocate (prdt_out  (npoin))
  allocate (hfx_out   (npoin))
  allocate (strs_out  (npoin))
  allocate(ust_out)
  allocate(sflux_out)
  allocate(LMO_out)
  ! allocate(T2_out)

end subroutine netCDF_allocate

subroutine netCDF_deallocate
  !
  !*** Close the file and deallocates
  !
  use def_master
  use def_netCDF
  use netCDF

  implicit none
  !  
  deallocate (time_out)
  deallocate (u_out)
  deallocate (v_out)
  deallocate (vel_out)
  deallocate (dir_out)
  deallocate (th_out)
  deallocate (uadv_out)
  deallocate (vadv_out)
  deallocate (thadv_out)
  deallocate (tke_out)
  deallocate (eps_out)
  deallocate (mut_out)
  deallocate (lm_out)
  deallocate (Rf_out)
  deallocate (prdm_out)
  deallocate (prdt_out)
  deallocate (hfx_out)
  deallocate (strs_out)
  deallocate(ust_out)
  deallocate(sflux_out)
  deallocate(LMO_out)
  ! deallocate(T2_out)

  return
end subroutine netCDF_deallocate
