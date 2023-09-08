subroutine tem_comput(itask,tloca,htota,cploc,temre)
  use def_parame
  implicit none
  integer(ip), intent(in) :: itask
  real(rp), intent(in)    :: tloca
  real(rp), intent(in)    :: htota
  real(rp), intent(in)    :: cploc(6,2)
  real(rp), intent(out)   :: temre
  real(rp)                :: hfute, htlow
  integer(ip)             :: trang 

  if (tloca < 1000.0_rp) then
    trang = 1_ip
  else
    trang = 2_ip
  end if

  select case(itask)

  case(1_ip)
    temre = (((cploc(5,trang)   * tloca + &
               cploc(4,trang) ) * tloca + &
               cploc(3,trang) ) * tloca + &
               cploc(2,trang) ) * tloca + &
               cploc(1,trang)

  case(2_ip)
    hfute = ((((cploc(5,trang) / 5.0_rp   * tloca + &
                cploc(4,trang) / 4.0_rp ) * tloca + &
                cploc(3,trang) / 3.0_rp ) * tloca + &
                cploc(2,trang) / 2.0_rp ) * tloca + &
                cploc(1,trang) )          * tloca + &
                cploc(6,trang)
    temre = hfute - htota

  end select

end subroutine tem_comput
