subroutine gabls3_ini
  !---------------------------------------------------------------
  !This subroutine applies for GABLS3 test case -> kfl_case.eq.3
  !
  ! 1st - Allocates specific case variables (mesoscale tendencies)
  ! 2nd - 
  !---------------------------------------------------------------
    
  use def_master
  implicit none

  real(8)      :: L,zp,velini,logs,tsfc,vel,zdif,shapes,eta,phi_h
  real(8)      :: z, kz, lm
  integer      :: i,j,thstar

  !Allocates specific GABLS3 variables
  allocate (tempe_adv(npoin))
  allocate (u_adv(npoin))
  allocate (v_adv(npoin))
  allocate (u_geo(npoin))
  allocate (v_geo(npoin))
  allocate (u_meso(npoin))
  allocate (v_meso(npoin))
  allocate (th_meso(npoin))

  ! Initial conditions interpolation for GABLS3
  i = 1
  do while (coord(ipoin).gt.height(i))
     i = i + 1
  end do
  ! Profiling from 1st Mesoscale Level down to Surface
  if ((i.eq.2.and.height(1).lt.1.0d-6).or.&    !height(1)=0
       (i.eq.1.and.height(1).gt.1.0d-6)) then  !height(1)>0
     ! Set first level
     if (height(1).lt.1.0d-6) then
        lev_ref = 2
        j    = 2
        zp   = height(2)
     else
        lev_ref = 1
        j    = 1
        zp   = height(1)
     end if
     ! Set Tempearture at surface (tsfc)
     velini = sqrt(uini(j)*uini(j)+vini(j)*vini(j))
     ustar  = (velini*kar)/log(1+(zp/rough))
     qw = hflux(1)*rhocp
     thstar = -qw/(ustar*rhocp)
     L = teref*ustar*ustar/(kar*gravi*thstar)
     logs = log(1+(2.0d0/rough))
     eta = 2.0d0/L
     if (kfl_temp.eq.1.or.kfl_temp.eq.3) then          !diagnosed from T2m
        !Stable              
        if (L.ge.0) then
           phi_h = mo_prand +mo_beta2*eta
           tsfc = t2(1)-((0.74d0*thstar/kar)*(logs+phi_h/mo_prand -1.0d0))             
        !Unstable
        else
           phi_h = mo_prand*(1.0 -mo_gamm2*eta)**(-0.5)
           tsfc = t2(1)-((0.74d0*thstar/kar)*(logs-2.0*log(0.5*(1.0+ mo_prand/phi_h))))                
        end if
     end if
     if (kfl_temp.eq.2) tsfc = tsk(2) !tsfc = tsk
     ! Velocity profiling down to Surface
     ! No Thermal Stability considered
     vel = (ustar/kar)*log(1.0d0+coord(ipoin)/rough)
     veloc_ini(ipoin,1) = vel*(uini(j)/velini)
     veloc_ini(ipoin,2) = vel*(vini(j)/velini)
     ! Temperature profiling down to Surface
     if (coord(ipoin).gt.2.0d0) then ! Above 2m
        zdif =  height(j)-2.0d0
        shapes = (height(j) - coord(ipoin))/zdif
        tempe_ini(ipoin)   = thini(j)*(1.0d0-shapes) + t2(1)*shapes
     else                            ! Below 2m          
        ! Stable
        logs = log(1+(coord(ipoin)/rough))
        if (L.ge.0.0d0) then
           phi_h = mo_prand +mo_beta2*eta
           tempe_ini(ipoin) = tsfc+((0.74d0*thstar/kar)*(logs+phi_h/mo_prand -1.0d0))              
        ! Unstable
        else
           phi_h = mo_prand*(1.0 -mo_gamm2*eta)**(-0.5)
           tempe_ini(ipoin) = tsfc+((0.74d0*thstar/kar)*(logs-2.0*log(0.5*(1.0+ mo_prand/phi_h))))
        end if
     end if

     ! At 0m height and Dirichlet Surface Boundary Condition
  else if  (ipoin.eq.1.and.coord(ipoin).eq.0.0d0.and.kfl_bouco_vel.eq.0) then
     ! Set first level
     if (height(1).lt.1.0d-6) then
        j    = 2
        zp   = height(2)
     else
        j    = 1
        zp   = height(1)
     end if
     ! Set Tempearture at surface (tsfc)
     if (kfl_temp.eq.1.or.kfl_temp.eq.3) then          !diagnosed from T2m
        velini = sqrt(uini(j)*uini(j)+vini(j)*vini(j))
        ustar  = (velini*kar)/log(1+(zp/rough))
        qw = hflux(1)*rhocp
        thstar = -qw/(ustar*rhocp)
        L = teref*ustar*ustar/(kar*gravi*thstar)
        logs = log(1+(2.0d0/rough))
        eta = 2/L
        !Stable              
        if (L.ge.0) then
           phi_h = mo_prand +mo_beta2*eta
           tsfc = t2(1)-((0.74d0*thstar/kar)*(logs+phi_h/mo_prand -1.0d0))              
        !Unstable
        else
           phi_h = mo_prand*(1.0 -mo_gamm2*eta)**(-0.5)
           tsfc = t2(1)-((0.74d0*thstar/kar)*(logs-2.0*log(0.5*(1.0+ mo_prand/phi_h))))               
        end if
     end if
     if (kfl_temp.eq.2) tsfc = tsk(2) !tsfc = tsk
     veloc_ini(1,1) = 0.0d0
     veloc_ini(1,2) = 0.0d0
     tempe_ini(1)   = tsfc

     ! Above 1st Mesoscale model -> Linear interpolation  
  else
     zdif =  height(i)-height(i-1)
     shapes = (height(i) - coord(ipoin))/zdif
     veloc_ini(ipoin,1) = uini(i)*(1.0d0-shapes) + uini(i-1)*shapes
     veloc_ini(ipoin,2) = vini(i)*(1.0d0-shapes) + vini(i-1)*shapes
     tempe_ini(ipoin)   = thini(i)*(1.0d0-shapes) + thini(i-1)*shapes
  end if
  
  z =  coord(ipoin)
  keyva_ini(ipoin)   = max(keyam,  ustar*ustar/sqrt(cmu) * (1.0d0 - z/1000.0 )**3.0) 
  kz =  kar*(z+rough)
  lm = kz/(1.0+ kz/l_max)
  epsil_ini(ipoin)   = max(epsam, ((cmu*keyva_ini(ipoin)*keyva_ini(ipoin))**0.75)/ lm)
  
end subroutine gabls3_ini



subroutine gabls2_ini
  !-------------------------------------------------------------------
  !This subroutine apply for GABLS2 test case -> kfl_case.eq. 1 and 2 
  !-------------------------------------------------------------------

  use def_master
  implicit none
  real(8)      :: z, kz, lm

  ctime = 16.0*3600.0 
  
  if (coord(ipoin).lt.200.0) then
     tempe_ini(ipoin) = 288.0 -2.0*coord(ipoin)/200.0
  else if (coord(ipoin).lt.850.0) then
     tempe_ini(ipoin) = 286.0
  else if (coord(ipoin).lt.900.0) then
     tempe_ini(ipoin) = 286.0 +2.0*(coord(ipoin)-850.0)/50.0
  else if (coord(ipoin).lt.1000.0) then
     tempe_ini(ipoin) = 288.0 +4.0*(coord(ipoin)-900.0)/100.0
  else if (coord(ipoin).lt.2000.0) then
     tempe_ini(ipoin) = 292.0 +8.0*(coord(ipoin)-1000.0)/1000.0
  else if (coord(ipoin).lt.3500.0) then
     tempe_ini(ipoin) = 300.0 +10.0*(coord(ipoin)-2000.0)/1500.0
  else if (coord(ipoin).lt.4000.0) then
     tempe_ini(ipoin) = 310.0 +2.0*(coord(ipoin)-3500.0)/500.0
  else
     tempe_ini(ipoin) = 312.0 +5.0*(coord(ipoin)-4000.0)/2000.0
  end if

  z =  coord(ipoin)
  keyva_ini(ipoin)   = max(keyam,  ustar*ustar/sqrt(cmu) * (1.0d0 - z/1000.0 )**3.0) 
  kz =  kar*(z+rough)
  lm = kz/(1.0+ kz/l_max)
  epsil_ini(ipoin)   = max(epsam, ((cmu*keyva_ini(ipoin)*keyva_ini(ipoin))**0.75)/ lm)     

end subroutine gabls2_ini



subroutine gabls3_begste
  !---------------------------------------------------------------
  ! This subroutine constructs BC's for GABLS3 case
  !---------------------------------------------------------------

  use def_master
  implicit none

  integer      :: i,j
  real(8)      :: dt,dz
  real(8)      :: ustar_wrf,thstar_wrf,L_wrf
  real(8)      :: h_shape, t_shape, before, after,logs,eta,phi_h
  real(8)      :: T_2,T_sfc,Th_lev,thstar,L,prandt
  real(8)      :: h_index(npoin)

  prandt = 0.74d0 

  !
  !Time shape function
  !
  i = 1
  do while (ctime.gt.seconds(i))
     i = i + 1
  end do
  !
  if (i.eq.1) then
     i = 2
     t_shape = 1.0d0
  else
     dt =  seconds(i)-seconds(i-1)
     t_shape = (seconds(i) - ctime)/dt
  end if

  !
  !Node height index for vertical interpolation
  !
  do ipoin = 1, npoin
     j = 1
     do while (coord(ipoin).gt.height(j))
        j = j + 1
     end do
     h_index(ipoin) = j
  end do
     
  !
  !Tewal diagnosed from T2m(WRF), qw(WRF) and ustar(WRF)
  !
  if (kfl_temp.eq.1) then
     T_2    =  t2(i)*(1.0d0-t_shape) + t2(i-1)*t_shape
     T_sfc  =  tsk(i)*(1.0d0-t_shape) + tsk(i-1)*t_shape
     ustar_wrf = ustini(i)*(1.0d0-t_shape) + ustini(i-1)*t_shape !WRF
     ! Calculates heat flux at wall
     qw = (hflux(i)*(1.0d0-t_shape) + hflux(i-1)*t_shape)*(rhocp)
     !
     thstar_wrf = -qw/(ustar_wrf*rhocp)
     L_wrf = teref*ustar_wrf*ustar_wrf/(kar*gravi*thstar_wrf)
     logs = log(1+(2.0d0/rough))
     eta = 2.0d0/L_wrf
     !
     !Stable
     if (L_wrf.ge.0.0d0) then
        phi_h = mo_prand +mo_beta2*eta
        tewal = T_2-((mo_prand*thstar_wrf/kar)*(logs+phi_h/mo_prand -1.0d0)) 
        write(13,'(f12.0,2x,5(f10.2,2x),e17.10,2x,4(f10.2,2x))') ctime,tewal,tempe(1,1),T_2,T_sfc,qw,L_wrf,ustar,ustar_wrf,1.0d0
    !Unstable
     else
        phi_h = mo_prand*(1.0 -mo_gamm2*eta)**(-0.5)
        tewal = T_2-((mo_prand*thstar_wrf/kar)*(logs-2.0*log(0.5*(1.0+ mo_prand/phi_h))))
        write(13,'(f12.0,2x,5(f10.2,2x),e17.10,2x,4(f10.2,2x))') ctime,tewal,tempe(1,1),T_2,T_sfc,qw,L_wrf,ustar,ustar_wrf,2.0d0
     end if

  !
  !Tewal diagnosed from T2m(WRF), qw(Wind1D) and ustar(Wind1D), (thstar(wind1d only))
  !
  else if (kfl_temp.eq.3) then
     ! Clipping for friction velocity
     if ((ustar.lt.0.05d0).or.(ustar2.lt.0.05d0)) then
        print *,'ustar,ustar2',ustar,ustar2
        ustar  = 0.05d0
        ustar2 = 0.05d0
        print *,'ustar and ustar2 set to 0.05'
     end if
     T_2    =  t2(i)*(1.0d0-t_shape) + t2(i-1)*t_shape
     T_sfc  =  tsk(i)*(1.0d0-t_shape) + tsk(i-1)*t_shape
     ! Calculates heat flux at wall
     if (istep.gt.1) then
        qw = -rhocp*ustar2*kar/(mo_prand*logft)*(tempe(1,1)-tewal)
     else
        qw = (hflux(i)*(1.0d0-t_shape) + hflux(i-1)*t_shape)*(rhocp)
     end if
     !
     thstar = -qw/(ustar2*rhocp)   !should I change ustar by ustar2?
     L = teref*ustar*ustar/(kar*gravi*thstar)
     logs = log(1+(2.0d0/rough))
     eta = 2.0d0/L
     !
     !Stable
     if (L.ge.0) then
        phi_h = mo_prand +mo_beta2*eta
        tewal = T_2-((mo_prand*thstar/kar)*(logs+phi_h/mo_prand -1.0d0)) 
        write(13,'(f12.0,2x,5(f10.2,2x),e17.10,2x,4(f10.2,2x))') ctime,tewal,tempe(1,1),T_2,T_sfc,qw,L,ustar,logft,1.0d0
     !Unstable
     else
        phi_h = mo_prand*(1.0 -mo_gamm2*eta)**(-0.5)
        tewal = T_2-((mo_prand*thstar/kar)*(logs-2.0*log(0.5*(1.0+ mo_prand/phi_h))))
        write(13,'(f12.0,2x,5(f10.2,2x),e17.10,2x,4(f10.2,2x))') ctime,tewal,tempe(1,1),T_2,T_sfc,qw,L,ustar,logft,2.0d0
     end if

  !   
  !Tewal equal to TSK(WRF)
  !
  else if (kfl_temp.eq.2) then
     Th_lev =  th(lev_ref,i)*(1.0d0-t_shape) + th(lev_ref,i-1)*t_shape
     T_2    =  t2(i)*(1.0d0-t_shape) + t2(i-1)*t_shape
     T_sfc  =  tsk(i)*(1.0d0-t_shape) + tsk(i-1)*t_shape
     tewal  = tsk(i)*(1.0d0-t_shape) + tsk(i-1)*t_shape
     write(13,'(f12.0,5(f8.2))') ctime,tewal,T_2,T_sfc,Th_lev,3.0d0
  else
     call runend('Transient temperature for GABLS3 not specified correctly')
  end if
  teref = tewal

  !
  ! Variable Top wind and Top temperature
  !
  !Top index
  j = h_index(npoin)
  !
  dz =  height(j)-height(j-1)
  h_shape = (height(j) - coord(npoin))/dz
  !U component Vertical and time interpolation
  before = u(j,i-1)*(1.0d0-h_shape) + u(j-1,i-1)*h_shape
  after  = u(j,i)*(1.0d0-h_shape) + u(j-1,i)*h_shape
  ugeos(1) = after*(1.0d0-t_shape) + before*t_shape
  !V component Vertical and time interpolation
  before = v(j,i-1)*(1.0d0-h_shape) + v(j-1,i-1)*h_shape
  after  = v(j,i)*(1.0d0-h_shape) + v(j-1,i)*h_shape
  ugeos(2) =  after*(1.0d0-t_shape) + before*t_shape
  !Variable top temperature
  before = th(j,i-1)*(1.0d0-h_shape) + th(j-1,i-1)*h_shape
  after  = th(j,i)*(1.0d0-h_shape) + th(j-1,i)*h_shape
  tetop =  after*(1.0d0-t_shape) + before*t_shape
  
  !
  ! Temperature Advection
  !
  if (kfl_thadv) then
     call interpola_3d(i,t_shape,h_index,thadv,tempe_adv)
  end if

  !
  ! Momentum Advections 
  !
  if (kfl_momadv) then
     call interpola_3d(i,t_shape,h_index,uadv,u_adv)
     call interpola_3d(i,t_shape,h_index,vadv,v_adv)
  end if

  !
  ! Geostrophic Velocity for Pressure gradient 
  !
  if (kfl_pressgr) then
     call interpola_3d(i,t_shape,h_index,ugeo,u_geo)
     call interpola_3d(i,t_shape,h_index,vgeo,v_geo)
  end if

  !
  ! Mesoscale Velocity and temperature for NUDGING
  !
  if (kfl_mom_nudging) then
     call interpola_3d(i,t_shape,h_index,u,u_meso)
     call interpola_3d(i,t_shape,h_index,v,v_meso)
  end if
  if (kfl_tem_nudging) then
     call interpola_3d(i,t_shape,h_index,th,th_meso)
  end if
  
end subroutine gabls3_begste


subroutine gabls2_begste(htime)
  ! This subroutine constructs BC's for GABLS2 case
  use def_master
  implicit none

  real(8),intent(in)   :: htime
  
  if (htime.lt.17.4d0) then
     tewal = -10.0d0 -25.0d0*cos(0.22d0*htime+0.2d0)
  else if (htime.lt.30d0) then
     tewal = -0.54d0*htime +15.2d0
  else if (htime.lt.41.9d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*htime+1.8d0)
  else if (htime.lt.53.3d0) then
     tewal = -0.37d0*htime +18.0d0
  else if (htime.lt.65.6d0) then
     tewal = -4.0d0 -25.0d0*cos(0.22d0*htime+2.5d0)
  else
     tewal = 4.4d0
  end if
  tewal = tewal + 273.15d0
  
end subroutine gabls2_begste


subroutine cycle_begste(htime)
  ! This subroutine constructs BC's for CYCLE based on GABLS2 case
  use def_master
  implicit none

  real(8),intent(in)   :: htime
  
  
  if (htime.lt.17.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime+24.0d0)+1.8d0)
  else if (htime.lt.29.8d0) then
     tewal = -0.37d0*(htime+24.0d0) +18.0d0
  else if (htime.lt.41.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*htime+1.8d0)
  else if (htime.lt.53.8d0) then
     tewal = -0.37d0*htime +18.0d0
  else if (htime.lt.65.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-24.0d0)+1.8d0)
  else if (htime.lt.77.8d0) then
     tewal = -0.37d0*(htime-24.0d0) +18.0d0
  else if (htime.lt.89.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-48.0d0)+1.8d0)
  else if (htime.lt.101.8d0) then
     tewal = -0.37d0*(htime-48.0d0) +18.0d0
  else if (htime.lt.113.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-72.0d0)+1.8d0)
  else if (htime.lt.125.8d0) then
     tewal = -0.37d0*(htime-72.0d0) +18.0d0
  else if (htime.lt.137.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-96.0d0)+1.8d0)
  else if (htime.lt.149.8d0) then
     tewal = -0.37d0*(htime-96.0d0) +18.0d0
  else if (htime.lt.161.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-120.0d0)+1.8d0)
  else if (htime.lt.173.8d0) then
     tewal = -0.37d0*(htime-120.0d0) +18.0d0
  else if (htime.lt.185.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-144.0d0)+1.8d0)
  else if (htime.lt.197.8d0) then
     tewal = -0.37d0*(htime-144.0d0) +18.0d0
  else if (htime.lt.209.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-168.0d0)+1.8d0)
  else if (htime.lt.221.8d0) then
     tewal = -0.37d0*(htime-168.0d0) +18.0d0
  else if (htime.lt.233.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-192.0d0)+1.8d0)
  else if (htime.lt.245.8d0) then
     tewal = -0.37d0*(htime-192.0d0)+18.0d0
  else if (htime.lt.257.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-216.0d0)+1.8d0)
  else if (htime.lt.269.8d0) then
     tewal = -0.37d0*(htime-216.0d0) +18.0d0
  else if (htime.lt.281.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-240.0d0)+1.8d0)
  else if (htime.lt.293.8d0) then
     tewal = -0.37d0*(htime-240.0d0) +18.0d0
  else if (htime.lt.305.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-264.0d0)+1.8d0)
  else if (htime.lt.317.8d0) then
     tewal = -0.37d0*(htime-264.0d0) +18.0d0
  else if (htime.lt.329.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-288.0d0)+1.8d0)
  else if (htime.lt.341.8d0) then
     tewal = -0.37d0*(htime-288.0d0) +18.0d0
  else if (htime.lt.353.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-312.0d0)+1.8d0)
  else if (htime.lt.365.8d0) then
     tewal = -0.37d0*(htime-312.0d0) +18.0d0
  else if (htime.lt.377.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-336.0d0)+1.8d0)
  else if (htime.lt.389.8d0) then
     tewal = -0.37d0*(htime-336.0d0) +18.0d0
  else if (htime.lt.401.935d0) then
     tewal = -7.0d0 -25.0d0*cos(0.21d0*(htime-360.0d0)+1.8d0)
  else
     tewal = 3.1386d0
  end if
  !
  tewal = tewal + 273.15d0

end subroutine cycle_begste

subroutine interpola_3d(i,t_shape,h_index,variable_in,variable_out)
  use def_master
  implicit none

  integer,intent(in)   :: i
  real(8),intent(in)   :: t_shape
  real(8),intent(in)   :: h_index(npoin)
  real(8),intent(in)   :: variable_in(nz,nt)
  real(8),intent(out)  :: variable_out(npoin)
  real(8)              :: h_shape, before, after
  integer              :: j
  real(8)              :: dz

  do ipoin = 1, npoin
     j = h_index(ipoin)
     if (j.eq.1) then
        before = variable_in(j,i-1)
        after  = variable_in(j,i)
        variable_out(ipoin) = after*(1.0d0-t_shape) + before*t_shape
     else
        dz =  height(j)-height(j-1)
        h_shape = (height(j) - coord(ipoin))/dz
        !Vertical and time interpolation
        before = variable_in(j,i-1)*(1.0d0-h_shape) + variable_in(j-1,i-1)*h_shape
        after  = variable_in(j,i)*(1.0d0-h_shape) + variable_in(j-1,i)*h_shape
        variable_out(ipoin) = after*(1.0d0-t_shape) + before*t_shape
     end if
  end do

end subroutine interpola_3d
