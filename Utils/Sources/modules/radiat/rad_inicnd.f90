subroutine rad_inicnd()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_inicnd
  ! NAME 
  !    rad_inicnd
  ! DESCRIPTION
  !    This routine sets up the initial condition for the radiation field
  !    according to the function KFL_INICO_RAD
  ! USED BY
  !    rad_iniunk
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_radiat
  use def_kermod
  use mod_ker_proper 
  implicit none
  integer(ip) :: ipoin,idime,icomp,ispec,dummi
  real(rp)    :: con,densi_rad(npoin),dummr(1)
 
  icomp=min(3_ip,ncomp_rad)

  select case(kfl_inico_rad)
     
  case(1)
     !
     ! Black body radiation
     !
     call ker_proper('DENSI','NPOIN',dummi,dummi,densi_rad)
     
     do ipoin=1,npoin
!!C        bvess_rad(ipoin,1)
!!$        con=0.0_rp
!!$        do ispec=1,nspec_rad
!!$           con = con+conce_rad(ipoin,ispec)*
!!$        enddo
        radav_rad(ipoin,icomp) = 4.0_rp * steph_rad * densi_rad(ipoin) *tempe_rad(ipoin)**4
     end do

  case(2)
     !
     ! Constant: Coefficients are given in reabcs
     !
     do ipoin=1,npoin
!!C        bvess_rad(ipoin,1)
        radav_rad(ipoin,icomp) = bvcoe_rad(1 ) !+kfl_fixno_rad(1,ipoin))  ! First coeff is for bulk, rest is for boundaries
     end do

  end select

end subroutine rad_inicnd
