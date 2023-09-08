subroutine rad_sgsope(itask,ielem,pgaus,gpsgs,gprad)     
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_sgsope
  ! NAME
  !   rad_sgsope
  ! DESCRIPTION
  !    Operate on the subgrid scale:
  !    ITASK = 1 ... Initialize SGS
  !          = 2 ... Add SGS to Radiation Intensity Average
  !          = 3 ... Update SGS
  ! OUTPUT
  ! USES
  ! USED BY
  !    rad_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master
  use def_radiat, only       :  kfl_sgsno_rad,resgs_rad,rasgs_rad
  implicit none
  integer(ip), intent(in)    :: itask,ielem,pgaus
  real(rp),    intent(inout) :: gprad(pgaus)
  real(rp),    intent(inout) :: gpsgs(pgaus,2)
  integer(ip)                :: igaus

  select case(itask)

  case(1)
     !
     ! Initialize SGS
     !
     if(kfl_sgsno_rad==1) then
        do igaus=1,pgaus
           gpsgs(igaus,1)=rasgs_rad(ielem)%a(igaus,1)
        end do        
     end if

  case(2)
     !
     ! Add SGS to radiation
     !
     if(kfl_sgsno_rad==1) then
        do igaus=1,pgaus
           gprad(igaus)=gprad(igaus)+gpsgs(igaus,1)
        end do
     end if

  case(3)
     !
     ! Compute SGS residual and update it
     !
     if(kfl_sgsno_rad==1) then
        do igaus=1,pgaus
           resgs_rad(1)=resgs_rad(1)+(gpsgs(igaus,1)-rasgs_rad(ielem)%a(igaus,1))**2
           resgs_rad(2)=resgs_rad(2)+gpsgs(igaus,1)**2
           rasgs_rad(ielem)%a(igaus,1)=gpsgs(igaus,1)
        end do
     end if

  end select

end subroutine rad_sgsope
