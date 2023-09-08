subroutine tem_calcEnthalpyBC()
!-----------------------------------------------------------------------
!****f* Temper/tem_calcEnthalpyBC
! NAME 
!    tem_endite
! DESCRIPTION
!    This routine re-computes the enthalpy after the properties were updated in chemic to 
!    allow obtaining correct temperature BC for isothermal cases.
! USES
!    tem_comput
! USED BY
!    tem_endite
!    tem_endste
!    tem_updbcs
!***
!-----------------------------------------------------------------------
  use def_kintyp,   only : ip,rp
  use def_master,   only : sphec,therm
  use def_domain,   only : elmar,npoin
  use def_temper,   only : kfl_fixno_tem,bvess_tem
  implicit none
  integer(ip) :: ipoin,ivalu
  real(rp)    :: dummr,tenew,cploc(6,2)

  do ipoin=1,npoin
    if(kfl_fixno_tem(1,ipoin)==1) then
      do ivalu = 1,6
        cploc(ivalu,1) = sphec(ipoin,ivalu,1)
        cploc(ivalu,2) = sphec(ipoin,ivalu,2)
      end do
      dummr = 0.0_rp
      call tem_comput(2_ip,bvess_tem(1,ipoin,2),dummr,cploc,tenew)
      bvess_tem(1,ipoin,1) = tenew
      therm(ipoin,1)       = bvess_tem(1,ipoin,1)
      therm(ipoin,2)       = bvess_tem(1,ipoin,1)
    end if
  end do


end subroutine tem_calcEnthalpyBC
