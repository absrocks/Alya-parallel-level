subroutine wav_elmsou(ndime,kfl_sourc_wav,sourc_wav,gpcod,cutim,gpsou)
  !------------------------------------------------------------------------
  !****f* Wavequ/wav_elmsou
  ! NAME 
  !    wav_elmsou
  ! DESCRIPTION
  !    Compute source terms
  ! USES
  ! USED BY
  !    wav_matrix
  !------------------------------------------------------------------------
  use      def_parame
  implicit none
  integer(ip), intent(in)  :: ndime,kfl_sourc_wav
  real(rp),    intent(in)  :: sourc_wav(*),gpcod(ndime),cutim
  real(rp),    intent(out) :: gpsou
  real(rp)                 :: fx,gt,t0,c1,c2,c3,piot2,piott

  select case(kfl_sourc_wav)

  case(-1)
     !
     ! Constant source term
     !
     gpsou = sourc_wav(1)

  case(1)
     !
     ! Ricker source (typical in geophysics):
     ! s(x,t) = f(x)g(t)
     !   f(x) = c1*exp(-c2*(x^2+y^2)
     !   g(t) = 2(pi/c3)^2*[2(pi/c3)^2*(t-c3)^2-1]*exp[-(pi/c3)^2*(t-c3)^2]
     !
     if(ndime==2) then
        fx    = sourc_wav(1)*exp(-sourc_wav(2)*(gpcod(1)*gpcod(1)+gpcod(2)*gpcod(2)))
     else
        fx    = sourc_wav(1)*exp(-sourc_wav(2)*(gpcod(1)*gpcod(1)+gpcod(2)*gpcod(2)&
             &                                 +gpcod(3)*gpcod(3)))        
     end if
     t0    = sourc_wav(3)
     piot2 = (pi/t0)*(pi/t0)
     piott = piot2*(cutim-t0)*(cutim-t0)
     gt    = 2.0_rp*piot2*(2.0_rp*piott-1.0_rp)*exp(-piott)
     gpsou = fx*gt

  end select

end subroutine wav_elmsou
