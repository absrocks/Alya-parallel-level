subroutine got_diffun(alpha,chale,diffu)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_diffun
  ! NAME 
  !    got_diffun
  ! DESCRIPTION
  !    This routine computes the diffusion DIFFU as:
  !
  !    DIFFU = u_max*alpha_max*C3*f(alpha) according to the function
  !    KFL_DIFUN_GOT:
  !
  !    1. Exp-log:   f = exp[ - log ( 1 + (alpha/(C1*alpha_max)) )^C2 ]
  !    2. Tanh:      f = 1 + tanh[ - (alpha/(C1*alpha_max))^C2 ] 
  !    3. Heaviside: f = 1 if alpha/alpha_max<C1
  !                    = 0 otherwise         
  !    or as follows 
  !           
  !    4. Element length: f = DIFFU = u_max*alpha_max*h*C1            
  ! USES
  ! USED BY
  !    got_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only  :  rp
  use def_master, only  :  ittim
  use def_gotita, only  :  diffu_got,almax_got,vemax_got,kfl_difun_got
  implicit none
  real(rp), intent(in)  :: alpha,chale(2)
  real(rp), intent(out) :: diffu
  real(rp)              :: funct

  select case(kfl_difun_got)

  case(1)
     !
     ! Exp-log function
     !
     funct=exp(-(log(alpha/(almax_got*diffu_got(1))+1.0_rp))**diffu_got(2))
     diffu=vemax_got*almax_got*diffu_got(3)*funct

  case(2)
     !
     ! Tanh function
     !
     funct=1.0+tanh(-(alpha/(almax_got*diffu_got(1)))**diffu_got(2))
     diffu=vemax_got*almax_got*diffu_got(3)*funct

  case(3)
     !
     ! Heaviside function
     !
     if(alpha/almax_got<diffu_got(1)) then
        funct=1.0_rp
     else
        funct=0.0_rp
     end if
     diffu=vemax_got*almax_got*diffu_got(3)*funct

  case(4)
     !
     ! Proportional to element length
     !
     diffu=vemax_got*almax_got*chale(1)*diffu_got(1)

  case(5)
     !
     ! Constant
     !
     diffu=vemax_got*almax_got*diffu_got(3)*diffu_got(1)

  end select

end subroutine got_diffun
