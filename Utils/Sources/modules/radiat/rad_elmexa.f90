subroutine rad_elmexa(&
     pgaus,gpcod,gpden,gpdif,gprea,gpgrd,&
     gpvel,gpsou,gprhs)
     
!-----------------------------------------------------------------------
!****f* Radiat/rad_elmexa
! NAME 
!    rad_elmexa
! DESCRIPTION
!    Compute RHS for exact solution
! USES
!    rad_exacso
! USED BY
!    rad_elmope
!***
!-----------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_domain, only    :  ndime
  use def_radiat, only    :  kfl_exacs_rad,expar_rad
  implicit none
  integer(ip), intent(in) :: pgaus
  real(rp),    intent(in) :: gpden(pgaus),gpdif(pgaus),gprea(pgaus)
  real(rp),    intent(in) :: gpgrd(ndime,pgaus),gpsou(pgaus)
  real(rp),    intent(in) :: gpvel(ndime,pgaus)
  real(rp),    intent(in) :: gpcod(ndime,pgaus)
  real(rp),    intent(in) :: gprhs(pgaus)
  real(rp)                :: dummr(3)
  integer(ip)             :: igaus

  if(kfl_exacs_rad/=0) then
     do igaus=1,pgaus
        call rad_exacso(&
             2_ip,gpcod(1,igaus),gpden(igaus),gpdif(igaus),&
             gprea(igaus),gpgrd(1,igaus),gpvel(1,igaus),gpsou(igaus),&
             dummr,dummr,gprhs(igaus))
     end do
  end if

end subroutine rad_elmexa
