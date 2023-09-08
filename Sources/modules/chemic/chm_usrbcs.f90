subroutine chm_usrbcs(iclas,itask)
  !-----------------------------------------------------------------------
  !****f* partis/chm_usrbcs
  ! NAME
  !    chm_usrbcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions
  ! OUTPUT 
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip), intent(in) :: iclas,itask
  integer(ip)             :: ipoin
  real(rp)                :: C0eq,Eform,Ceq,ovkT,T

  if( iclas < 1 .or. iclas > nclas_chm ) call runend('CHM_USRBCS: WRONG CLASS')

  select case(itask)

  case(1_ip)

     bvess_chm(iclas,278)=1e5

  case(2_ip)

     call chm_usrtem(T)
     ovkT   = 1.0_rp/(boltz_chm*T) 
     C0eq   = equil_chm(1,iclas)                 ! C0eq
     Eform  = equil_chm(2,iclas)                 ! E_form: formation energy
     Ceq    = C0eq*exp(-Eform*ovkT)              ! Ceq = C0eq * exp(-E_form/kT)
     do ipoin=1,npoin
        bvess_chm(iclas,ipoin)=Ceq
     end do

  case(3_ip)

     do ipoin = 1,npoin
        bvess_chm(iclas,ipoin) = 2.0_rp * coord(1,ipoin) + 3.0_rp
     end do

  end select

end subroutine chm_usrbcs
