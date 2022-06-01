subroutine write_restart
  !------------------------------------------------
  ! This routine writes restart files
  !------------------------------------------------

  use def_master
  implicit none
  character*40 :: file
  real(8)      :: gpnut,tstar
  integer      :: ipoin
  !
  !** Write Single value varibales
  !
  gpnut =  cmu*keyva(1,1)*keyva(1,1)/epsil(1,1)
  hflx0 =  rhocp*gpnut*(tempe(2,1)- tempe(1,1))/(coord(2)-coord(1))/sigte
 
  file='Single_value.rst'
  open(lun_rst1, FILE=file, status='unknown' )
  write(lun_rst1,'(a1,9a19)')  '#', '1:istep', '2:time', '3:nodes', '4:ustar', '5:heatfl', &
       '6:tewal', '7:maxle', '8:fcori', '9:ustar2'
  write (lun_rst1,'(i17,2x,e22.15,2x,i17,2x,6(e22.15,2x))') istep, ctime, npoin, ustar, hflx0, tewal, &
       lenmy, fcori, ustar2
  close (lun_rst1)

  !
  !** Write profiles variables
  !
  file='Profiles.rst'
  open(lun_rst2, FILE=file, status='unknown' )
  write(lun_rst2,'(a,6a17)')  '#', '1:coord', '2:U', '3:V', '4:key', '5:eps', '6:temp'
  do ipoin = 1, npoin
     write (lun_rst2, '(6(e22.15, x))') coord(ipoin),veloc(ipoin,1,1),veloc(ipoin,2,1), &
          keyva(ipoin,1), epsil(ipoin,1), tempe(ipoin,1)
  end do
  close (lun_rst2)

end subroutine write_restart



subroutine read_restart
  !----------------------------------------------------
  ! This routine interpolates restart variables to node
  !----------------------------------------------------

  use def_master
  implicit none

  character*40       :: file, cvoid
  integer            :: iz,i,j,step, ipoin, ielem
  logical            :: mesh_coincide = .false., prev_thermal = .true.
  real(8)            :: diff,L,vel,velini,logs,phi_h,eta,shapes,zdif,thstar,ustar2_rst, time_rst
  real(8)            :: tnume, tdeno,  sqkey, gpcod, dz, z
  !
  !*** Reads restart files
  !
  ! Reads Single value variables
  file = 'Single_value.rst'  ! plotprof_hh.mm
  open(90, FILE=trim(file),status='unknown')
  read(90,*) cvoid
  !1:istep 2:time 3:nodes 4:ustar 5:heatfl 6:tewal 7:maxle 8:fcori '9:ustar2''
  read(90,'(i17,2x,e22.15,2x,i17,2x,6(e22.15,2x))') step,time_rst, nz_rst, ust_rst, hflux_rst, tewal_rst, lmax_rst, fcori_rst,ustar2_rst 
  close(90)
  ! Checks Coriolis Force
  if (abs(fcori_rst-fcori).gt.1.0d-9) write(*,*) 'RESTART: Different Coriolis force',fcori_rst,fcori

  !
  ! Overwrite Single value variables
  
  
  ustar  = ust_rst
  ustar2 = ustar2_rst
  logft =  log(1.0+ dwall/rough)
  !
  ! Allocates Porfiles variables
  allocate(z_rst(nz_rst))
  allocate(u_rst(nz_rst))
  allocate(v_rst(nz_rst))
  allocate(key_rst(nz_rst))
  allocate(eps_rst(nz_rst))
  allocate(temp_rst(nz_rst))

  !
  ! Reads Profiles variables
  file = 'Profiles.rst'  ! plotprof_hh.mm
  open(90, FILE=trim(file),status='unknown')
  read(90,*) cvoid
  do iz=1,nz_rst
     read(90,'(6(e22.15, x))') z_rst(iz),u_rst(iz),v_rst(iz),key_rst(iz),eps_rst(iz),temp_rst(iz)
  end do
  close(90)
  if (maxval(temp_rst).lt.1.0d-4) prev_thermal = .false.
  if (.not.prev_thermal) then
     ctime = 0.0d0
     istep = 0
     hflx0 = 0.0d0
  else
     istep  = step
     ctime  = time_rst
     tewal  = tewal_rst
     teref  = tewal
     lenmy  = lmax_rst
     hflx0  = hflux_rst
     
  end if
  ! Checks mesh coincidence
  if (npoin.eq.nz_rst) then
     diff = 0.0d0
     do i=1,npoin
        diff = diff + abs(coord(i)-z_rst(i))
     end do
     if (diff.lt.1.d-03) mesh_coincide = .true.
  end if

  !
  !*** Set Restart variables as Initial Conditions 
  !
  if (mesh_coincide) then !dump data
     write (*,*) 'RESTART: mesh coincide'
     veloc_ini(:,1) = u_rst(:) 
     veloc_ini(:,2) = v_rst(:)
     keyva_ini(:)   = key_rst(:) 
     epsil_ini(:)   = eps_rst(:)
     if (prev_thermal)  tempe_ini(:)   = temp_rst(:)   
     !
  else !Linear Interpolation
     write (*,*) 'RESTART: mesh does not coincide'
     j = 1 ! Restart levels counter
     do i=1,npoin
        !
        !Coordinate below first restart node
        if ((coord(i)-z_rst(1)).lt.-1.0d-6) then   

           ! Velocity profiling down to Surface
           ! No Thermal Stability considered
           velini = sqrt(u_rst(1)*u_rst(1)+v_rst(1)*v_rst(1))
           vel = (ustar/kar)*log(1.0d0+coord(i)/rough)    
           veloc_ini(i,1) = vel*(u_rst(1)/velini)
           veloc_ini(i,2) = vel*(v_rst(1)/velini)

           ! Temperature profiling down to Surface using MO
           ! Thermal Stability considered
           thstar =  hflx0/(rhocp*ustar)
           L = teref*ustar*ustar/(kar*gravi*thstar)
           eta = z_rst(1)/L
           logs = log(1+(coord(i)/rough))
           if (L.ge.0.0d0) then ! Stable
              phi_h = mo_prand +mo_beta2*eta
              if (prev_thermal)  tempe_ini(i) = tewal+((0.74d0*thstar/kar)*(logs+phi_h/mo_prand -1.0d0))              
           else ! Unstable
              phi_h = mo_prand*(1.0 -mo_gamm2*eta)**(-0.5)
              if (prev_thermal) tempe_ini(i) = tewal+((0.74d0*thstar/kar)*(logs-2.0*log(0.5*(1.0+ mo_prand/phi_h))))
           end if

           ! TKE and Dissipation equal to first restart node
           keyva_ini(i) = key_rst(1)
           epsil_ini(i) = eps_rst(1)

        else ! Above first restart node

           ! Looks for the reference node           
           do while ((coord(i)-z_rst(min(j,nz_rst))).ge.1.0d-6)
              j = j + 1
              ! Node above top height restart file
              if (j.gt.nz_rst) then 
                 write(*,*) 'RESTART: Node above top height restart file coord=',coord(i)
                 veloc_ini(i,1) = u_rst(nz_rst) 
                 veloc_ini(i,2) = v_rst(nz_rst)
                 keyva_ini(i)   = key_rst(nz_rst) 
                 epsil_ini(i)   = eps_rst(nz_rst)
                 if (prev_thermal) tempe_ini(i)   = temp_rst(nz_rst)
                 exit
              end if
           end do
           
           ! Interpolation
           if (j.eq.1) j = 2

           if (j.le.nz_rst) then
              zdif =  z_rst(j)-z_rst(j-1)
              shapes = (z_rst(j) - coord(i))/zdif
              veloc_ini(i,1) = u_rst(j)*(1.0d0-shapes) + u_rst(j-1)*shapes
              veloc_ini(i,2) = v_rst(j)*(1.0d0-shapes) + v_rst(j-1)*shapes
              keyva_ini(i)   = key_rst(j)*(1.0d0-shapes) + key_rst(j-1)*shapes
              epsil_ini(i)   = eps_rst(j)*(1.0d0-shapes) + eps_rst(j-1)*shapes
              if (prev_thermal) tempe_ini(i)   = temp_rst(j)*(1.0d0-shapes) + temp_rst(j-1)*shapes
           end if
        end if
     end do
  end if

  !
  !*** Allocates specific GABLS3 problem type variables
  !
  if (kfl_case.eq.3) then
     allocate (tempe_adv(npoin))
     allocate (u_adv(npoin))
     allocate (v_adv(npoin))
     allocate (u_geo(npoin))
     allocate (v_geo(npoin))
     allocate (u_meso(npoin))
     allocate (v_meso(npoin))
     allocate (th_meso(npoin))
  end if

  ! set initial temp profile. 
  if (.not.prev_thermal) then
  !   tempe_ini = tewal
     teref = tewal
     tnume = 0.0d0
     tdeno = 0.0d0
     do ielem =1, nelem
        ipoin = ielem
        jpoin = ielem + 1
        dz = coord(jpoin)-coord(ipoin)
        sqkey = 0.50d0*(sqrt(keyva_ini(ipoin))+sqrt(keyva_ini(jpoin)))
        gpcod = coord(ipoin) +0.50*dz
        tnume = tnume + gpcod*sqkey*dz
        tdeno = tdeno + sqkey*dz
     end do
     lenmy = 0.075*tnume/tdeno
     do ipoin =1, npoin
        z = coord(ipoin)
        if (z.lt.ztemin) then
           tempe_ini(ipoin) =  tewal +  gradbo*z
        else
           tempe_ini(ipoin) =  tewal +  gradbo*ztemin + gradto*(z-ztemin)
        end if
     end do
  end if
  
  
  return
end subroutine read_restart
