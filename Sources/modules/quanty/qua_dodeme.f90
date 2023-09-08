subroutine qua_dodeme(order)
  !------------------------------------------------------------------------
  !****f* Quanty/qua_dodeme
  ! NAME 
  !    qua_dodeme
  ! DESCRIPTION
  !    This subroutine is the bridge between Quanty module and 
  !    Dodeme service
  ! USES
  ! USED BY
  !    qua_matrix
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_quanty
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: order

  if(kfl_servi(7)/=0) then

     select case(order)

     case(1)
        !
        ! Ask Dodeme if natural transmission conditions
        ! must be computed. If this is the case NINTE_DOD=1
        !
        ivari_dod=30+ID_TEMPE
        ninte_dod=0
        call Dodeme(9_ip)
        if(ninte_dod/=0) then
           !
           ! Compute natural transmission conditions
           !
!           bemol_dod=bemol_tem
!           call tem_diffus(diffu_dod) 
!           if(kfl_advec_tem==1) then
!              veloc_dod => veloc(:,:,1)
!           else
!              call Dodeme(11_ip)
!              call tem_velfun(kfl_advec_tem,ndime,npoin,coord,veloc_dod)
!           end if
           call Dodeme(10_ip) 
!           if(kfl_advec_tem==1) nullify(veloc_dod)
           !
           ! Assemble natural transmission conditions
           !
           call Dodeme(8_ip)
        end if
        !
        ! Ask Dodeme if integral b.c.
        ! must be computed. If this is the case NINTE_DOD=1
        !
        ivari_dod=ID_TEMPE
        ninte_dod=0
        call Dodeme(9_ip)
        if(ninte_dod/=0) then
           !
           ! Compute essential projection transmission conditions
           !
!           kfl_fixno_dod => kfl_fixno_qua(1,:)
!           call Dodeme(10_ip) 
!           nullify(kfl_fixno_dod)
        end if
        !
        ! Assemble nodal transmission contribution
        !
        call Dodeme(8_ip)

     case(2)

     end select

  end if

end subroutine qua_dodeme
