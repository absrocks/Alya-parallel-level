subroutine ibm_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_updunk
  ! NAME 
  !    ibm_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates 
  ! USED BY
  !    ibm_endste (itask=5)
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use mod_kdtree

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iimbo,idime,ipoib
  real(rp)                :: xcoor(3),rotco(3)
  real(rp),    pointer    :: F(:,:),a(:,:),v(:,:),x(:,:)          ! Linear motion
  real(rp),    pointer    :: T(:,:),z(:,:),w(:,:),s(:,:),quate(:,:) ! Angular motion

  select case (itask)

  case(1_ip)
     
     !-------------------------------------------------------------------
     !
     ! Initial accelerations
     ! 
     !-------------------------------------------------------------------

      do iimbo = 1,nimbo      
        a =>  imbou(iimbo) % accel
        v =>  imbou(iimbo) % velol
        x =>  imbou(iimbo) % posil
        z =>  imbou(iimbo) % accea
        w =>  imbou(iimbo) % veloa
        s =>  imbou(iimbo) % posia
        z =>  imbou(iimbo) % accea
        if( kfl_linib_ibm /= 0 ) then
           do idime = 1,ndime
              a(idime,2) = a(idime,1)
              v(idime,2) = v(idime,1)
              x(idime,2) = x(idime,1)
           end do
        end if
        if( kfl_rotib_ibm /= 0 ) then
           do idime = 1,3
              w(idime,2) = w(idime,1)
              s(idime,2) = s(idime,1)
              z(idime,2) = z(idime,1)
           end do
        end if
     end do
  case(5_ip)

     !-------------------------------------------------------------------
     !
     ! End of time step
     ! 
     !-------------------------------------------------------------------

     do iimbo = 1,nimbo
        F =>  imbou(iimbo) % force
        a =>  imbou(iimbo) % accel
        v =>  imbou(iimbo) % velol
        x =>  imbou(iimbo) % posil

        T =>  imbou(iimbo) % torqu
        z =>  imbou(iimbo) % accea
        w =>  imbou(iimbo) % veloa
        s =>  imbou(iimbo) % posia
        quate => imbou(iimbo) % quate
        !
        ! New time step: n <= n + 1
        !
        quate(1,2) = quate(1,1)
        do idime = 1,3
           F(idime,2)   = F(idime,1) 
           a(idime,2)   = a(idime,1) 
           v(idime,2)   = v(idime,1) 
           x(idime,2)   = x(idime,1)

           T(idime,2)   = T(idime,1) 
           z(idime,2)   = z(idime,1) 
           w(idime,2)   = w(idime,1) 
           s(idime,2)   = s(idime,1)                 
           quate(idime+1,2) = quate(idime+1,1)
        end do
     end do


  end select

end subroutine ibm_updunk

