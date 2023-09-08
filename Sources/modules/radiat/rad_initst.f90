subroutine rad_initst()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_initst
  ! NAME 
  !    rad_initst
  ! DESCRIPTION
  !    This routine initializes test cases temperature and parameters
  ! USED BY
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  implicit none
  integer(ip)             :: ipoin,ispec
  real(rp)                :: radius

  if( INOTMASTER ) then
     
     select case (kfl_atest_rad)
     case(1_ip)   ! Case 1: Cylinder with hot gas inside
        aniso_rad(1)=0.0_rp
        scatt_rad(1)=0.0_rp
        absor_rad(1)=expar_rad(1) 
        do ipoin=1,npoin
           do ispec=1,nspec_rad
              conce_rad(ipoin,ispec) = 1.0_rp
           enddo
           if (coord(3,ipoin) <= 0.0_rp .or. coord(3,ipoin) >= 2.0_rp .or. coord(1,ipoin)**2+coord(2,ipoin)**2 >= 1.0_rp) then
              tempe_rad(ipoin) = 0.0_rp   ! Zero temp walls :)
           else
              tempe_rad(ipoin)=1000.0_rp  !Hot core
           endif
        enddo
     case(2_ip)   ! Rectangular plane with periodic boundaries in x
        aniso_rad(1)=0.0_rp
        scatt_rad(1)=0.0_rp
        absor_rad(1)=expar_rad(1) 
        do ipoin=1,npoin
           do ispec=1,nspec_rad
              conce_rad(ipoin,ispec) = 1.0_rp
           enddo
           tempe_rad(ipoin)=(100.0_rp**4+(200.0_rp**4-100.0_rp**4)*coord(2,ipoin))**0.25
        enddo
     case(3_ip)   ! Two concentric spheres
        aniso_rad(1)=0.0_rp
        scatt_rad(1)=0.0_rp
        absor_rad(1)=expar_rad(1) 
        do ipoin=1,npoin
           conce_rad(ipoin,1) = 1.0_rp
           radius = sqrt(coord(1,ipoin)**2+coord(2,ipoin)**2+coord(3,ipoin)**2)
           if (radius.le.0.25_rp) then !inner shell
              tempe_rad(ipoin)=10.0_rp
           else if(radius.ge.0.5_rp) then !Top wall
              tempe_rad(ipoin)=50.0_rp
           else
              tempe_rad(ipoin)= 10.0_rp+(50.0_rp-10.0_rp)*(radius-0.25_rp)/(0.5_rp-0.25_rp)! Linear Bulk temp
           endif
        enddo
     end select

  endif

end subroutine rad_initst
