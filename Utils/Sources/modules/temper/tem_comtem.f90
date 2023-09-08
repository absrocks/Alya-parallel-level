subroutine tem_comtem
  !-----------------------------------------------------------------------
  !****f* Temper/solve_tem
  ! NAME 
  !    solve_tem
  ! DESCRIPTION
  !    This routine computes the temperature from the enthalpy using the 
  !    polynomial coefficients for the specific heat
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_ker_proper
  implicit none
  integer(ip)             :: ipoin,dummi,nte,xt, ivalu
  real(rp)                :: dummy(ndime,ndime),teloc(2),delta
  real(rp)                :: cploc(6,2),deltaDash,ent_loc,acval

  !
  ! total enthalpy equation, source term implicit in cp coefficients
  !
  do ipoin=1,npoin
    !
    ! Check if adiabatic calculation
    !
    if (kfl_adiab_tem /= 0) then           ! Enthalpy is given by the mixture fraction
       acval   = min(1.0_rp,max(0.0_rp,(conce(ipoin,3,1))))
       ent_loc = acval * cfi_hmax_tem + (1-acval) * cfi_hmin_tem
    else
       ent_loc = therm(ipoin,1)            ! Enthalpy is recomputed each time-step
    endif
    teloc(1) = tempe(ipoin,1)

    nte = 0
    do ivalu = 1,6
      cploc(ivalu,1) = sphec(ipoin,ivalu,1)
      cploc(ivalu,2) = sphec(ipoin,ivalu,2)
    end do

    do while(nte < 50)
      call tem_comput(2_ip,teloc(1),ent_loc,cploc,delta)
      !
      nte = nte + 1
      if (abs(delta) < 0.00001_rp) then
        !
        ! update tempe and cp and exit loop
        !
        tempe(ipoin,1) = min(max(teloc(1),200.0_rp),3000.0_rp)
        call tem_comput(1_ip,tempe(ipoin,1),ent_loc,cploc,sphek(ipoin,1))
        !
        exit
      else
        !
        ! residual still too high, compute new temperature for next loop
        !
        call tem_comput(1_ip,teloc(1),ent_loc,cploc,deltaDash)
        !
        teloc(2) = teloc(1)
        teloc(1) = teloc(2) - (delta/deltaDash)
      end if
    end do
  end do

end subroutine tem_comtem
