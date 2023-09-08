subroutine chm_reamet()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_reamet
  ! NAME 
  !    chm_reamet
  ! DESCRIPTION
  !    This routine reads meteo file (different formats are possible)
  !    Meteo data comes from (mesoscale) meteo models
  ! USES
  !    
  ! USED BY
  !    chemic
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_kintyp
  use def_master
  use def_domain
  use def_chemic
  use mod_memchk  
  use def_elmtyp
  implicit none
!
  logical(lg)           :: go_on
  character(20)         :: wline
  integer(ip)           :: ipoin,idime,iline
  integer(ip)           :: ixmet,iymet,izmet,ivmet
  real(rp)              :: time1,time2,dummr
!
  integer(ip)           :: ndims_met(6)
  integer(ip)           :: ix_met,iy_met,iz_met(4)
  real(rp)              :: shapf(8) 
  real(rp), allocatable :: var2d_met(:,:,:), var3d_met(:,:,:,:) 
                         
!
!  Debug option 
!
  if ( kfl_meteo_chm == 1 ) then      
      !
      if( INOTMASTER ) then
      !
      !  Sequential or slave
      !
        do ipoin = 1,npoin
           veloc_chm(    1,ipoin) = 10.0_rp
           veloc_chm(    2,ipoin) = 0.0_rp
           veloc_chm(ndime,ipoin) = 0.0_rp
           tempe_chm(      ipoin) = 1.0_rp
           densi_chm(      ipoin) = 1.0_rp
        end do
      end if
      !
      tmete_chm = 1e10
      return
  end if
!
!  Other options
!
!  1. Sequential or master. Read the meteo file (in the required format) for the current time slab.
!      
!     NOTE: It is assumed that meteo variables are on a structured mesh regular on XY plane. However, values for
!           z (vertical) can vary.
!
!     Meteo variables are stored in arrays to minimize subsequent communication during broadcast calls.
!     INTEGER
!         ndims_met(1) = nx_met    
!         ndims_met(2) = ny_met 
!         ndims_met(3) = nz_met 
!         ndims_met(4) = 2d_met          Number of 2D met variables read
!         ndims_met(5) = 3d_met          Number of 3D met variables read
!         ndims_met(6) = INT(tmete_chm)  
!     REAL2D
!         var2d_met(nx_met,ny_met,2d_met) where
!         var2d_met(:,:,1)     x_met
!         var2d_met(:,:,2)     y_met   
!     REAL3D
!         var3d_met(nx_met,ny_met,nz_met,3d_met) where
!         var3d_met(:,:,:,1)   z_met
!         var3d_met(:,:,:,2)   u_met
!         var3d_met(:,:,:,3)   v_met
!         var3d_met(:,:,:,4)   w_met
!         var3d_met(:,:,:,5)   t_met
!         var3d_met(:,:,:,6)   r_met
!
  if( INOTSLAVE ) then
     !
     if ( kfl_meteo_chm == 2 ) then 
       !
       !  ASCII file. File format is
       !        time1 time2
       !        nx_met,ny_met,nz_met
       !          x_met y_met                              nx_met*ny_met        values              
       !          ...
       !          ...
       !          z_met u_met v_met w_met t_met r_met      nx_met*ny_met*nz_met values
       !          ...
       !          ...
       !
       ndims_met(4) = 2   ! 2D
       ndims_met(5) = 6   ! 3D
       !
       go_on = .true.
       iline = 0
       do while(go_on)
          iline = iline+1
          read(lun_remet_chm,*,err=101) time1,time2      
          iline = iline+1
          read(lun_remet_chm,*,err=101) ndims_met(1),ndims_met(2),ndims_met(3)
          if( oltim >= time1 .and. oltim <= time2 ) then
              !
              !  Time interval found
              !
              go_on = .false.
              allocate(var2d_met(ndims_met(1),ndims_met(2),             ndims_met(4)))
              allocate(var3d_met(ndims_met(1),ndims_met(2),ndims_met(3),ndims_met(5)))
              do iymet = 1,ndims_met(2)
              do ixmet = 1,ndims_met(1)
                 iline = iline+1
                 read(lun_remet_chm,*,err=101) (var2d_met(ixmet,iymet,ivmet),ivmet=1,ndims_met(4))
              end do
              end do
              do izmet = 1,ndims_met(3)
              do iymet = 1,ndims_met(2)
              do ixmet = 1,ndims_met(1)
                 iline = iline+1
                 read(lun_remet_chm,*,err=101) (var3d_met(ixmet,iymet,izmet,ivmet),ivmet=1,ndims_met(5))
             end do
             end do
             end do
         else
             !
             !  Time interval not found
             !   
             do ipoin = 1, ndims_met(1)*ndims_met(2)
                iline = iline+1
                read(lun_remet_chm,*,err=101) dummr
             end do
             do ipoin = 1, ndims_met(1)*ndims_met(2)*ndims_met(3)
                iline = iline+1
                read(lun_remet_chm,*,err=101) dummr
             end do
         end if
     end do
     ndims_met(6)  = INT(time2)
!
!    
    else if ( kfl_meteo_chm == 3 ) then 
     !
     !  netCDF    
     !
      call runend('chm_reamet: netcdf not implemented yet')
     !
    end if
!
!
!
 end if
!  
!
!  2. In parallel flow, master sends information to slaves
!
!
  if( IPARALL ) then
    !
    call parari('BCT',0_ip,SIZE(ndims_met),ndims_met) 
    !
    if(ISLAVE) then
       allocate(var2d_met(ndims_met(1),ndims_met(2),             ndims_met(4)))
       allocate(var3d_met(ndims_met(1),ndims_met(2),ndims_met(3),ndims_met(5)))
    end if
    !
    call pararr('BCT',0_ip,ndims_met(1)*ndims_met(2)*ndims_met(4)             ,var2d_met) 
    call pararr('BCT',0_ip,ndims_met(1)*ndims_met(2)*ndims_met(3)*ndims_met(5),var3d_met) 
   !
  end if
!
!
!  3. Sequential or slave. Linear interpolation of meteo variables at each ipoin
!
!
  if( INOTMASTER ) then
   !
    do ipoin = 1,npoin
    !
    !  Determinates interpolation parameters for each point 
    !
       call chm_intmet(coord(1,ipoin),coord(2,ipoin),coord(3,ipoin),var2d_met(1,1,1),var2d_met(1,1,2),var3d_met(1,1,1,1), &
                      ndims_met(1),ndims_met(2),ndims_met(3),ix_met,iy_met,iz_met,shapf)
    !
    !  Interpolates
    ! 
       do idime = 1,ndime
          veloc_chm(idime,ipoin) = shapf(1) * var3d_met(ix_met  ,iy_met  ,iz_met(1)  ,idime+1) + &
                                   shapf(2) * var3d_met(ix_met+1,iy_met  ,iz_met(2)  ,idime+1) + &
                                   shapf(3) * var3d_met(ix_met+1,iy_met+1,iz_met(3)  ,idime+1) + &
                                   shapf(4) * var3d_met(ix_met  ,iy_met+1,iz_met(4)  ,idime+1) + &
                                   shapf(5) * var3d_met(ix_met  ,iy_met  ,iz_met(1)+1,idime+1) + &
                                   shapf(6) * var3d_met(ix_met+1,iy_met  ,iz_met(2)+1,idime+1) + &
                                   shapf(7) * var3d_met(ix_met+1,iy_met+1,iz_met(3)+1,idime+1) + &
                                   shapf(8) * var3d_met(ix_met  ,iy_met+1,iz_met(4)+1,idime+1)
       end do
       tempe_chm(ipoin) = shapf(1) * var3d_met(ix_met  ,iy_met  ,iz_met(1)  ,5) + &
                          shapf(2) * var3d_met(ix_met+1,iy_met  ,iz_met(2)  ,5) + &
                          shapf(3) * var3d_met(ix_met+1,iy_met+1,iz_met(3)  ,5) + &
                          shapf(4) * var3d_met(ix_met  ,iy_met+1,iz_met(4)  ,5) + &
                          shapf(5) * var3d_met(ix_met  ,iy_met  ,iz_met(1)+1,5) + &
                          shapf(6) * var3d_met(ix_met+1,iy_met  ,iz_met(2)+1,5) + &
                          shapf(7) * var3d_met(ix_met+1,iy_met+1,iz_met(3)+1,5) + &
                          shapf(8) * var3d_met(ix_met  ,iy_met+1,iz_met(4)+1,5)
       densi_chm(ipoin) = shapf(1) * var3d_met(ix_met  ,iy_met  ,iz_met(1)  ,6) + &
                          shapf(2) * var3d_met(ix_met+1,iy_met  ,iz_met(2)  ,6) + &
                          shapf(3) * var3d_met(ix_met+1,iy_met+1,iz_met(3)  ,6) + &
                          shapf(4) * var3d_met(ix_met  ,iy_met+1,iz_met(4)  ,6) + &
                          shapf(5) * var3d_met(ix_met  ,iy_met  ,iz_met(1)+1,6) + &
                          shapf(6) * var3d_met(ix_met+1,iy_met  ,iz_met(2)+1,6) + &
                          shapf(7) * var3d_met(ix_met+1,iy_met+1,iz_met(3)+1,6) + &
                          shapf(8) * var3d_met(ix_met  ,iy_met+1,iz_met(4)+1,6)
    end do

   ! 
   !  Exchange between slaves for potential sources at the border elements
   !
   !   call pararr('SLX',NPOIN_TYPE,ndime*npoin,veloc_chm)
   !   call pararr('SLX',NPOIN_TYPE,      npoin,tempe_chm)
   !   call pararr('SLX',NPOIN_TYPE,      npoin,densi_chm)
!
  end if
!
!
!  4. Release memeory
!
    tmete_chm = 1.0_rp*ndims_met(6) 
    deallocate(var2d_met)
    deallocate(var3d_met)
!
    return

101 continue
    wline = intost(iline)
    call runend('ERROR WHILE READING CHEMIC METEO FILE AT LINE '//trim(wline))
!
end subroutine chm_reamet
! 
!
!
subroutine chm_intmet(x,y,z,x_met,y_met,z_met,nx_met,ny_met,nz_met,ix_met,iy_met,iz_met,shapf)
 !-----------------------------------------------------------------------
  !****f* chemic/chm_intmet
  ! NAME 
  !    chm_intmet
  ! DESCRIPTION
  !    This routine computes the interpolation factors for the point of
  !    coordinates (x,y,z) 
  !    OUTPUT
  !       ix_met,iy_met,iz_met(4)
  !       shapf(8)
  ! USES
  !    
  ! USED BY
  !    chm_reamet
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_kintyp
  use def_elmtyp
  use mod_elmgeo, only  :  elmgeo_natural_coordinates 
  implicit none
!
  integer(ip) :: nx_met,ny_met,nz_met,ix_met,iy_met,iz_met(4)
  real   (rp) :: x_met(nx_met,ny_met),y_met(nx_met,ny_met),z_met(nx_met,ny_met,nz_met)
  real   (rp) :: x,y,z,shapf(8)
!
  logical(lg) :: go_on
  integer(ip) :: ip_met,ii,jj,ifoun
  real   (rp) :: elcod(2,4),shapff(4),deriv(2,4),coglo(3),coloc(3)
  real   (rp) :: s,t,w(4),sm,tm,sp,tp
!
!  First, search on XY plane (it is assumed that the meteo data mesh is layered) for ix_met and iy_met
!
   coglo(1) = x
   coglo(2) = y
   go_on = .true.
   ip_met = 0
   do while(go_on)
      ip_met = ip_met + 1
      iy_met = (ip_met-1)/nx_met + 1
      ix_met = ip_met - (iy_met-1)*nx_met
!
      ix_met = min(ix_met,nx_met-1)
      iy_met = min(iy_met,ny_met-1)
!
      elcod(1,1) = x_met(ix_met  ,iy_met  )
      elcod(1,2) = x_met(ix_met+1,iy_met  )
      elcod(1,3) = x_met(ix_met+1,iy_met+1)
      elcod(1,4) = x_met(ix_met  ,iy_met+1)
      elcod(2,1) = y_met(ix_met  ,iy_met  )
      elcod(2,2) = y_met(ix_met+1,iy_met  )
      elcod(2,3) = y_met(ix_met+1,iy_met+1)
      elcod(2,4) = y_met(ix_met  ,iy_met+1)
      !
      !  If the element is found reutrns ifoun=1 and returns shapf,deriv,coloc
      !
      call elmgeo_natural_coordinates(      &
           2_ip,QUA04,4_ip,elcod,shapff,   &
           deriv,coglo,coloc,ifoun)
      !
      if(ifoun == 1) then
         go_on = .false.
         s = coloc(1)
         t = coloc(2)
      end if
!
      if( (ip_met == nx_met*ny_met).and.(go_on) ) then 
         call runend(' chm_intmet : interpolation not found on xy plane')
      end if  
!
 end do
!
! Second, search on Z plane for iz_met(4)
!
   do ip_met = 1,4
      iz_met(ip_met) = 0
      if(ip_met.eq.1) then
         ii = ix_met
         jj = iy_met
      else if(ip_met.eq.2) then
         ii = ix_met + 1
         jj = iy_met
      else if(ip_met.eq.3) then
         ii = ix_met + 1
         jj = iy_met + 1
      else if(ip_met.eq.4) then
         ii = ix_met 
         jj = iy_met + 1
      end if
      !
      if(z.le.z_met(ii,jj,1)) then
      !
      !  bottom
      !
         iz_met(ip_met) = 1
         w(ip_met)      = 0.0            
      else if(z.ge.z_met(ii,jj,nz_met)) then
      !
      !  top
      !
         iz_met(ip_met) = nz_met - 1
         w(ip_met)      = 1.0    
      else
      !
      !  middle
      !
         go_on = .true.
         do while(go_on)
            iz_met(ip_met) = iz_met(ip_met) + 1
            if(z.ge.z_met(ii,jj,iz_met(ip_met)) .and. z.le.z_met(ii,jj,iz_met(ip_met)+1)) then
                 go_on = .false. 
                 w(ip_met) = (z-z_met(ii,jj,iz_met(ip_met)))/ & 
                             (z_met(ii,jj,iz_met(ip_met)+1)-z_met(ii,jj,iz_met(ip_met)))
             end if
         end do
       end if
   end do
!
!  Computes interpolation shapf functions
!
   s = max(s,-1.0_rp)       ! parameter s in (-1,1)
   s = min(s, 1.0_rp)
   t = max(t,-1.0_rp)       ! parameter s in (-1,1)
   t = min(t, 1.0_rp)
   w = max(w, 0.0_rp)       ! parameter w in (0,1)
   w = min(w, 1.0_rp)
!
   sm = 0.5_rp*(1.0_rp-s)
   tm = 0.5_rp*(1.0_rp-t)
   sp = 0.5_rp*(1.0_rp+s)
   tp = 0.5_rp*(1.0_rp+t)
!      
   shapf(1) = (1.0_rp-w(1))*sm*tm
   shapf(2) = (1.0_rp-w(2))*sp*tm
   shapf(3) = (1.0_rp-w(3))*sp*tp
   shapf(4) = (1.0_rp-w(4))*sm*tp
   shapf(5) = (       w(1))*sm*tm
   shapf(6) = (       w(2))*sp*tm
   shapf(7) = (       w(3))*sp*tp
   shapf(8) = (       w(4))*sm*tp
!
  return
  end subroutine chm_intmet
