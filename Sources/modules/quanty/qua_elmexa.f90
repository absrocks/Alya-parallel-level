subroutine qua_elmexa(&
     pgaus,gpcod,gpden,gpdif,gprea,gpgrd,&
     gpvel,gpsou,gprhs)
     
!-----------------------------------------------------------------------
!****f* Quanty/qua_elmexa
! NAME 
!    qua_elmexa
! DESCRIPTION
!    Compute RHS for exact solution
! USES
!    qua_exacso
! USED BY
!    qua_elmope
!***
!-----------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_domain, only    :  ndime
  use def_quanty

  implicit none
  integer(ip), intent(in) :: pgaus
  real(rp),    intent(in) :: gpden(pgaus),gpdif(pgaus),gprea(pgaus)
  real(rp),    intent(in) :: gpgrd(ndime,pgaus),gpsou(pgaus)
  real(rp),    intent(in) :: gpvel(ndime,pgaus)
  real(rp),    intent(in) :: gpcod(ndime,pgaus)
  real(rp),    intent(in) :: gprhs(pgaus)
  real(rp)                :: dummr(3)
  integer(ip)             :: igaus

  if(kfl_exacs_qua/=0) then
     do igaus=1,pgaus
        call qua_exacso(&
             2_ip,gpcod(1,igaus),gpden(igaus),gpdif(igaus),&
             gprea(igaus),gpgrd(1,igaus),gpvel(1,igaus),gpsou(igaus),&
             dummr,dummr,gprhs(igaus))
     end do
  end if

end subroutine qua_elmexa
