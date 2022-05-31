subroutine ibm_eulers(dt)
  !-----------------------------------------------------------------------
  !****f* ibm_eulers/ibm_eulers
  ! NAME
  !    ibm_eulers
  ! DESCRIPTION
  !    This routines solves the Euler and Newton equations for rigid bodies
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_immbou
  use mod_kdtree
  use mod_messages, only : livinf

  implicit none
  real(rp), intent(in) :: dt
  integer(ip)          :: iimbo,idime
  real(rp)             :: onvma,auxir
  real(rp)             :: Faver(3),Taver(3)

  real(rp), pointer    :: F(:,:),a(:,:),v(:,:)        ! Linear motion
  real(rp), pointer    :: T(:,:),z(:,:),IT(:),R(:,:)  ! Angular motion

  if( ittim /= 0 ) then
     if( kfl_colli_ibm == 0 ) then
        call livinf(165_ip,'EULER,',0_ip) 
     else
        call livinf(165_ip,'EULER',0_ip) 
        call livinf(164_ip,' ',0_ip) 
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! First time step
  !
  !----------------------------------------------------------------------
  if( ittim == 1 ) then      
     do iimbo = 1,nimbo        
        onvma = 1.0_rp / imbou(iimbo)%massa       
        F =>  imbou(iimbo) % force
        T =>  imbou(iimbo) % torqu
        if( kfl_linib_ibm /= 0 ) then
           do idime = 1,ndime     
              F(idime,2) = F(idime,1) 
           end do
        end if
        if( kfl_rotib_ibm /= 0 ) then
           do idime = 1,3
              T(idime,2) = T(idime,1)
           end do
        end if
     end do
  end if
  !----------------------------------------------------------------------
  !
  ! Solver Euler equations
  ! 
  !----------------------------------------------------------------------

  do iimbo = 1,nimbo

     onvma = 1.0_rp / imbou(iimbo)%massa

     F =>  imbou(iimbo) % force
     a =>  imbou(iimbo) % accel
     v =>  imbou(iimbo) % velol


     T =>  imbou(iimbo) % torqu
     z =>  imbou(iimbo) % accea
    
     Faver(3) = 0.0_rp
     Taver(3) = 0.0_rp
     if (kfl_nforc_ibm == 1) then
        do idime = 1,3
           Faver(idime) = F(idime,1)
           Taver(idime) = T(idime,1)
        end do
     elseif (kfl_nforc_ibm == 2) then
        do idime = 1,3
           Faver(idime) = 0.5_rp*(F(idime,1) + F(idime,2))
           Taver(idime) = 0.5_rp*(T(idime,1) + T(idime,2))
        end do
     end if
     !----------------------------------------------------------------------
     !     
     ! Compute linear quantities
     !
     !----------------------------------------------------------------------
     
     
     if( kfl_linib_ibm /= 0 ) then
        do idime = 1,ndime
           if ( kfl_ralei_ibm == 0 ) then          
              a(idime,1) =              &
                   xline_ibm(idime) * Faver(idime) * onvma
           else

              if ( ittim <= nstra_ibm ) then 
                 auxir = 1.0_rp
              else if ( ittim >= nenra_ibm ) then 
                 auxir = 0.0_rp
              else
                 auxir = dble ( nenra_ibm - ittim ) / dble ( nenra_ibm - nstra_ibm )
              end if

              a(idime,1) =              &
                   xline_ibm(idime) * ( Faver(idime) * onvma - &
                   ( auxir * ralei_ibm * v(idime,2) ) )
           end if

           if( ittim == 1 ) a(idime,2) = a(idime,1)
        end do
     end if

     !----------------------------------------------------------------------
     !     
     ! Compute angular quantities
     !
     !----------------------------------------------------------------------
     R => imbou(iimbo) % rotac
     IT => imbou(iimbo)%momin
     if( kfl_rotib_ibm /= 0 ) then
        if ( ndime == 2 ) then
           !
           ! Compute angular quantities in 2D
           !
           z(3,1) = xrota_ibm(3) * Taver(3) / IT(1)        
           if( ittim == 1 ) z(3,2) = z(3,1)
        else if ( ndime == 3 ) then
           call ibm_fixpo1(dt,iimbo,0.000001_rp,100_ip)
        end if
     end if

  end do

end subroutine ibm_eulers
