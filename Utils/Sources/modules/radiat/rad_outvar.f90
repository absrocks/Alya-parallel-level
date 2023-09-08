subroutine rad_outvar(ivari)
  !------------------------------------------------------------------------
  !****f* Radiat/rad_output
  ! NAME 
  !    rad_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    rad_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  use mod_gradie

  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: ibopo,ipoin,idime
  real(rp)                :: dummr,rutim
  real(rp)                :: rad_gamma

  if( ivari == 0 ) return
  !
  ! Define postprocess variable
  !
  rutim = cutim

  select case (ivari)  

  case(1_ip)
     !
     ! Radiation Average
     !
     gesca => radav_rad(:,1) 

  case(2_ip)
     !
     ! Radiation Heat source term
     !
     gesca => radso(:,1) 

  case(3_ip)
     !
     ! Error w/r exact solution
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call rad_exaerr(2_ip)
     end if

  case(4_ip)
     !
     ! Residual
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = radav_rad(ipoin,1)-radav_rad(ipoin,2)
        end do
     end if
     rutim = real(ittot_rad)

   case(5_ip)
      !
      ! Heat flux vector
      !
      if( INOTMASTER ) then
         call memgen(zero,ndime,npoin)
         call gradie(radav_rad(:,1),gevec)
         do ipoin = 1,npoin
            do idime = 1,ndime
               gevec(idime,ipoin) = - rad_gamma(absor_rad(1), scatt_rad(1), aniso_rad(1)) * gevec(idime,ipoin)
            enddo
         end do
      end if
  case(6_ip)
     !
     ! Internal temperature used by RADIAT
     !
     gesca => tempe_rad
 
  end select

  call outvar(&
       ivari,&
       ittim,rutim,postp(1) % wopos(1,ivari))

end subroutine rad_outvar
