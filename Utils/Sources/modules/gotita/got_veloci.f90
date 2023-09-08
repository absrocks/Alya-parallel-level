subroutine got_veloci(ipoin,xcoor,xvelo)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_veloci
  ! NAME 
  !    got_veloci
  ! DESCRIPTION
  !    Returns the velocity
  ! USES
  ! USED BY
  !    got_turnon
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_master, only     :  veloc
  use def_gotita, only     :  kfl_velfu_got,veloc_got
  implicit none
  integer(ip), intent(in)  :: ipoin
  real(rp),    intent(in)  :: xcoor(ndime)
  real(rp),    intent(out) :: xvelo(ndime)
  integer(ip)              :: idime

  select case(kfl_velfu_got)

  case(-1)
     !
     ! Uses VELOC
     !
     do idime=1,ndime
        xvelo(idime) = veloc_got(idime,ipoin)
     end do

  case(0)
     !
     ! Uses VELOC
     !
     do idime=1,ndime
        xvelo(idime) = veloc(idime,ipoin,1)
     end do

  case(1)
     !
     ! u=(1,1)
     !
     do idime=1,ndime
        xvelo(idime)=1.0_rp
     end do

  case(2)
     !
     ! u=(x,-y)
     !
     xvelo(1)= xcoor(1)
     xvelo(2)=-xcoor(2)

  case(3)
     !
     ! u=(1,0)
     !
     xvelo(1)=1.0_rp
     xvelo(2)=0.0_rp

  case default
     !
     ! Error
     !
     call runend('UNKNOWN AIR VELOCITY FUNCTION')

  end select

end subroutine got_veloci
