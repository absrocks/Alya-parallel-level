subroutine ibm_movbox(iimbo,iposi,tstep,bobox)

  use def_kintyp, only    : ip,rp
  use def_domain, only    : ndime
  use def_master, only    : imbou,dtime
  use def_immbou
  implicit none

  integer(ip), intent(in)  :: iimbo
  integer(ip), intent(in)  :: iposi
  real(rp),    intent(in)  :: tstep
  real(rp),    intent(out) :: bobox(3,2)
  integer(ip)              :: idime
  real(rp)                 :: xi(3),si(3),smaxi,tempo
  !
  ! Calcute the linear acceleration and the angular and linear displacement in the next time step. We suppose a linear acceleration
  ! We use the Taylor Series for the approximation
  !  
  do idime = 1,3
     xi(idime) = tstep*imbou(iimbo)%velol(idime,iposi) + 0.5_rp*tstep*tstep*imbou(iimbo)%accel(idime,3) + &
          (1.0_rp/(6.0_rp*dtime))*tstep*tstep*tstep*(imbou(iimbo)%accel(idime,1) - imbou(iimbo)%accel(idime,2))
     si(idime) = tstep*imbou(iimbo)%veloa(idime,iposi) + 0.5_rp*tstep*tstep*imbou(iimbo)%accea(idime,3) +&
          (1.0_rp/(6.0_rp*dtime))*tstep*tstep*tstep*(imbou(iimbo)%accea(idime,1) - imbou(iimbo)%accea(idime,2))           
  end do
  
  smaxi = 0.0_rp
  do idime=1,3
     smaxi = smaxi + si(idime)*si(idime)
  end do
  smaxi = sqrt(smaxi)
  smaxi = smaxi*imbou(iimbo)%maxdi
  
  do idime = 1,3
     tempo = max( 0.0_rp, imbou(iimbo)%maxdi - 0.5_rp*(imbou(iimbo)%bobox(idime,2) - imbou(iimbo)%bobox(idime,1)) )
     tempo = min( tempo , smaxi)
     bobox(idime,1) = min( imbou(iimbo)%bobox(idime,1), imbou(iimbo)%bobox(idime,1) + xi(idime) - tempo )
     bobox(idime,2) = max( imbou(iimbo)%bobox(idime,2), imbou(iimbo)%bobox(idime,2) + xi(idime) + tempo )
  end do
  
end subroutine ibm_movbox






